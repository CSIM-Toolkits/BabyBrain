#include "itkImageFileWriter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryContourImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkNeighborhoodIterator.h"
#include "itkConstantBoundaryCondition.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkMedianImageFilter.h"
#include "itkMaskImageFilter.h"

#include "itkPluginUtilities.h"

#include "BrainVolumeRefinementCLP.h"

#include <ctime>

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{

template <typename TPixel>
int DoIt( int argc, char * argv[], TPixel )
{
    PARSE_ARGS;

    std::clock_t begin = clock();

    typedef TPixel InputPixelType;
    typedef unsigned short LabelPixelType;

    const unsigned int Dimension = 3;

    typedef itk::Image<InputPixelType,  Dimension>    InputImageType;
    typedef itk::Image<LabelPixelType, Dimension>     LabelImageType;

    typedef itk::ImageFileReader<InputImageType>  ImageReaderType;

    typename ImageReaderType::Pointer inputReader = ImageReaderType::New();

    inputReader->SetFileName( inputVolume.c_str() );
    inputReader->Update();

    typedef itk::LabelStatisticsImageFilter<InputImageType, LabelImageType>             LabelStatisticsImageFilter;
//    typedef itk::BinaryThresholdImageFilter<InputImageType, InputImageType>             BinaryThresholdInputType;
    typedef itk::BinaryThresholdImageFilter<InputImageType, LabelImageType>             BinaryThresholdType;
    typedef itk::BinaryContourImageFilter<LabelImageType, InputImageType>               BinaryContourType;
    typedef itk::SubtractImageFilter<InputImageType, InputImageType>                    SubtractLabelType;
    typedef itk::MedianImageFilter<InputImageType, InputImageType>                      MedianFilterType;
    typedef itk::ConnectedComponentImageFilter<LabelImageType, LabelImageType>          ConnectedLabelType;
    typedef itk::RelabelComponentImageFilter<LabelImageType, LabelImageType>            RelabelerType;
    typedef itk::BinaryThresholdImageFilter<LabelImageType,LabelImageType>              BinaryLabelMap;
    typedef itk::MaskImageFilter<InputImageType, LabelImageType>                        MaskFilterType;

    typename MedianFilterType::Pointer smoothInput = MedianFilterType::New();
    typename MedianFilterType::RadiusType mRadius;
    if (doInitialFillHoles) {
        //The user can choose if a initial hole filling procedure is passed into the input image (based on a global median filtering)
        std::cout<<"Correction of the initial brain mask requested..."<<std::endl;
        smoothInput->SetInput(inputReader->GetOutput());
        mRadius[0] = medianRadius[0];
        mRadius[1] = medianRadius[1];
        mRadius[2] = medianRadius[2];
        smoothInput->SetRadius(mRadius);
        smoothInput->Update();
    }
    std::cout<<"Selection mode: "<<selectionMode<<std::endl;

    // Copying information and pixel values from the input image
    typename InputImageType::Pointer updatedImage = InputImageType::New();
    updatedImage->CopyInformation(inputReader->GetOutput());
    updatedImage->SetRegions(inputReader->GetOutput()->GetRequestedRegion());
    updatedImage->Allocate();
    updatedImage->FillBuffer(0);
    itk::ImageRegionIterator<InputImageType>    inIt(((doInitialFillHoles)?smoothInput->GetOutput():inputReader->GetOutput()), inputReader->GetOutput()->GetRequestedRegion());
    itk::ImageRegionIterator<InputImageType>    upIt(updatedImage, updatedImage->GetRequestedRegion());
    upIt.GoToBegin();
    inIt.GoToBegin();
    while (!upIt.IsAtEnd()) {
        upIt.Set(inIt.Get());
        ++upIt;
        ++inIt;
    }

    //Creating auxiliary image to assist the iterative procedure
    typename InputImageType::Pointer auxImage = InputImageType::New();
    auxImage->CopyInformation(inputReader->GetOutput());
    auxImage->SetRegions(inputReader->GetOutput()->GetRequestedRegion());
    auxImage->Allocate();
    auxImage->FillBuffer(0);
    itk::ImageRegionIterator<InputImageType>    auxIt(auxImage, auxImage->GetRequestedRegion());
    auxIt.GoToBegin();
    inIt.GoToBegin();
    while (!auxIt.IsAtEnd()) {
        auxIt.Set(inIt.Get());
        ++auxIt;
        ++inIt;
    }

    //If user's want it, the new brain mask is returned.
    typename LabelImageType::Pointer newBrainMask = LabelImageType::New();
    newBrainMask->CopyInformation(inputReader->GetOutput());
    newBrainMask->SetRegions(inputReader->GetOutput()->GetRequestedRegion());
    newBrainMask->Allocate();
    newBrainMask->FillBuffer(0);

    //Auxiliary skeleton to infer the searching path for each iteration
    typename InputImageType::Pointer previousSkeleton = InputImageType::New();
    previousSkeleton->CopyInformation(inputReader->GetOutput());
    previousSkeleton->SetRegions(inputReader->GetOutput()->GetRequestedRegion());
    previousSkeleton->Allocate();
    previousSkeleton->FillBuffer(0);

    //Auxiliary image to infer the actual cutting region of the image
    typename InputImageType::Pointer cuttingArea = InputImageType::New();
    cuttingArea->CopyInformation(inputReader->GetOutput());
    cuttingArea->SetRegions(inputReader->GetOutput()->GetRequestedRegion());
    cuttingArea->Allocate();
    cuttingArea->FillBuffer(0);



    std::cout<<"*********************"<<std::endl;
    std::cout<<"BrainVolumeRefinement: Initiated"<<std::endl;
    std::cout<<"Iteration: ";
    bool stoppedByConvergence = false;
    for (int n = 0; n < numberOfIterations; ++n) {
        typename BinaryThresholdType::Pointer estimateMask = BinaryThresholdType::New();
        typename BinaryContourType::Pointer maskContour = BinaryContourType::New();
        typename SubtractLabelType::Pointer actualContour = SubtractLabelType::New();
        typename BinaryThresholdType::Pointer newMask = BinaryThresholdType::New();
        typename ConnectedLabelType::Pointer connLabel = ConnectedLabelType::New();
        typename RelabelerType::Pointer relabel = RelabelerType::New();
        typename BinaryLabelMap::Pointer connectedMask = BinaryLabelMap::New();
        typename MaskFilterType::Pointer maskInput = MaskFilterType::New();

        typename LabelImageType::SizeType radius; //The radius of the neighborhood used as the searching path

        //Calculating the brain volume of iteration i
        double totalVolumeI;
        itk::ImageRegionIterator<InputImageType>    itI(updatedImage, updatedImage->GetRequestedRegion());
        itI.GoToBegin();
        while (!itI.IsAtEnd()) {
            if (itI.Get()!=static_cast<InputPixelType>(0)) {
                totalVolumeI++;
            }
            ++itI;
        }

        //Making the image mask of the input volume
        std::cout<<(n+1)<<"..."<<std::flush;
        estimateMask->SetInput(updatedImage);
        estimateMask->SetLowerThreshold(1.0);
        estimateMask->SetUpperThreshold(itk::NumericTraits<InputPixelType>::max());
        estimateMask->SetInsideValue( 1 );
        estimateMask->Update();

        //************************************************************************************************************************
        //OPTIONAL: Dependes on the selection Mode
        typename LabelStatisticsImageFilter::Pointer globalMaskStatistics = LabelStatisticsImageFilter::New();
        globalMaskStatistics->SetLabelInput( estimateMask->GetOutput() );
        globalMaskStatistics->SetInput(updatedImage);
        globalMaskStatistics->Update();

        InputPixelType globalMean = static_cast<InputPixelType>(0);
        globalMean = (globalMaskStatistics->GetMean( 1 ))/static_cast<InputPixelType>(2); //TODO Como melhorar essa suposicao de dividir por 2
        //************************************************************************************************************************
        //Create volume contour from the input mask
        maskContour->SetInput(estimateMask->GetOutput());
        maskContour->SetForegroundValue( 1 );
        maskContour->Update();

        //Calculating the actual skeleton for iteration i and its major region, this is used to infer the searching region area which is used to clean the wrong voxels in the outside part of image.
        actualContour->SetInput1(maskContour->GetOutput());
        actualContour->SetInput2(previousSkeleton);
        actualContour->Update();

        //Searching region. This is useful to restrict the actual region to a minor part of the image, which decrease the number of points that the iteratior will visit.
        LabelImageType::SizeType windowSize;
        LabelImageType::IndexType initIndex;
        itk::ImageRegionConstIteratorWithIndex<InputImageType> modifyRegionIt(actualContour->GetOutput(),actualContour->GetOutput()->GetRequestedRegion());
        int minX=itk::NumericTraits<int>::max(), maxX=0;
        int minY=itk::NumericTraits<int>::max(), maxY=0;
        int minZ=itk::NumericTraits<int>::max(), maxZ=0;
        modifyRegionIt.GoToBegin();
        while (!modifyRegionIt.IsAtEnd()) {
            if (modifyRegionIt.Get()>static_cast<InputPixelType>(0)) {
                //Finding minima
                //X index
                if (modifyRegionIt.GetIndex()[0] < minX) {
                    minX = modifyRegionIt.GetIndex()[0];
                }
                //Y index
                if (modifyRegionIt.GetIndex()[1] < minY) {
                    minY = modifyRegionIt.GetIndex()[1];
                }
                //Z index
                if (modifyRegionIt.GetIndex()[2] < minZ) {
                    minZ = modifyRegionIt.GetIndex()[2];
                }

                //Finding maxima
                //X index
                if (modifyRegionIt.GetIndex()[0] > maxX) {
                    maxX = modifyRegionIt.GetIndex()[0];
                }
                //Y index
                if (modifyRegionIt.GetIndex()[1] > maxY) {
                    maxY = modifyRegionIt.GetIndex()[1];
                }
                //Z index
                if (modifyRegionIt.GetIndex()[2] > maxZ) {
                    maxZ = modifyRegionIt.GetIndex()[2];
                }
            }
            ++modifyRegionIt;
        }
        initIndex[0] = minX;
        initIndex[1] = minY;
        initIndex[2] = minZ;
        windowSize[0] = (maxX - minX);
        windowSize[1] = (maxY - minY);
        windowSize[2] = (maxZ - minZ);

        LabelImageType::RegionType searchingRegion;
        if (updatedImage->GetRequestedRegion().GetSize()[0]>windowSize[0] ||
                updatedImage->GetRequestedRegion().GetSize()[1]>windowSize[1] ||
                updatedImage->GetRequestedRegion().GetSize()[2]>windowSize[2]) {
            //If the estimated region is bigger than the original image, then do not change it (stopped by using the convergence criteria).
            searchingRegion.SetIndex(initIndex);
            searchingRegion.SetSize(windowSize);
        }else{
            stoppedByConvergence=true;
            break;
        }

        std::cout<<"Searching region: "<<searchingRegion<<std::endl;
        radius[0] = neighborRadius[0]; // TODO Fazer uma estimativa do tamanho do raio a partir do spacing (1,1,1 pode ser 3,3,3, mas 0.2,0.2,1.0 poderia ser 15,15,3 para manter a mesma proporcao de pixel sendo recrutados... deixar flag bool no XML para isto)
        radius[1] = neighborRadius[1];
        radius[2] = neighborRadius[2];
        itk::ImageRegionIterator<InputImageType> contourIt(actualContour->GetOutput(),searchingRegion);
        itk::NeighborhoodIterator<InputImageType, itk::ConstantBoundaryCondition<InputImageType> > imageIt(radius, updatedImage,searchingRegion);
        itk::NeighborhoodIterator<InputImageType, itk::ConstantBoundaryCondition<InputImageType> > updateIt(radius, auxImage,searchingRegion);

        itk::NeighborhoodIterator<InputImageType, itk::ConstantBoundaryCondition<InputImageType> > setForRemoveIt(radius, cuttingArea,searchingRegion);

        /* Running over the both image and mask contour in order to find the closest voxel in the brain tissue.
         * The local mean intensity, calculated over a defined neighborhood, is used as voting for removel criteria
         * which creates marks on another image (cuttingArea) that will be used for further local voxel removal. */
        contourIt.GoToBegin();
        imageIt.GoToBegin();
        updateIt.GoToBegin();

        setForRemoveIt.GoToBegin();
        unsigned int N = (neighborRadius[0]*2 + 1)*(neighborRadius[1]*2 + 1)*(neighborRadius[2]*2 + 1);
        InputPixelType wIntensity;
        while (!imageIt.IsAtEnd()) {
            if (contourIt.Get()>static_cast<LabelPixelType>(0)) {
                //Calculating local statistics into the neighborhood

                if (selectionMode=="Local") {
                    wIntensity=0.0;
                    for (unsigned int p = 0; p < N; p++) {
                        //Weighted average to local intensity
                        wIntensity += imageIt.GetPixel(p);
                    }
                    wIntensity/=N;

                    //Cutting out voxels that does not belongs to the brain based on the mean gray level intensity
                    for (unsigned int p = 0; p < N; p++) {
                        if (imageIt.GetPixel(p)!=static_cast<InputPixelType>(0)) {
                            if (imageIt.GetPixel(p) < (wIntensity)){
                                //Voting voxels that helps us to infer the outside area of the brain
                                setForRemoveIt.SetPixel(p,1);
                            }
                        }
                    }
                }else if (selectionMode=="Global") {
                    //Cutting out voxels that does not belongs to the brain based on the mean gray level intensity
                    for (unsigned int p = 0; p < N; p++) {
                        if (imageIt.GetPixel(p)!=static_cast<InputPixelType>(0)) {
                            if (imageIt.GetPixel(p) < (globalMean)){
                                //Voting voxels that helps us to infer the outside area of the brain
                                setForRemoveIt.SetPixel(p,1);
                            }
                        }
                    }
                }else if (selectionMode=="Manual") {
                    //Cutting out voxels that does not belongs to the brain based on the mean gray level intensity
                    for (unsigned int p = 0; p < N; p++) {
                        if (imageIt.GetPixel(p)!=static_cast<InputPixelType>(0)) {
                            if (imageIt.GetPixel(p) < static_cast<InputPixelType>(manualThreshold)){
                                //Voting voxels that helps us to infer the outside area of the brain
                                setForRemoveIt.SetPixel(p,1);
                            }
                        }
                    }
                }

            }

            ++contourIt;
            ++imageIt;
            ++updateIt;
            ++setForRemoveIt;
        }

        // Running over the image frontier again in order to ...
        //        radius[0] = neighborRadius[0]; // TODO Colocar uma variavel para o tamanho da jaenal de corte
        //        radius[1] = neighborRadius[1];
        //        radius[2] = neighborRadius[2];
        itk::NeighborhoodIterator<InputImageType, itk::ConstantBoundaryCondition<InputImageType> > imageCuttingIt(radius, updatedImage,searchingRegion);
        itk::NeighborhoodIterator<InputImageType, itk::ConstantBoundaryCondition<InputImageType> > setForRemoveCuttingIt(radius, cuttingArea,searchingRegion);

        contourIt.GoToBegin();
        imageCuttingIt.GoToBegin();
        setForRemoveCuttingIt.GoToBegin();

        //        unsigned int totalCuttingN = (radius[0]*2 + 1)*(radius[1]*2 + 1)*(radius[2]*2 + 1);
        while (!imageCuttingIt.IsAtEnd()) {
            //            int countImgVoxels=0, countRemovingVoxels=0;
            if (setForRemoveCuttingIt.GetCenterPixel()>static_cast<LabelPixelType>(0)) {
                //                //Cutting out voxels that does not belongs ...
                //                for (unsigned int p = 0; p < totalCuttingN; p++) {
                //                    if (imageCuttingIt.GetPixel(p)!=static_cast<InputPixelType>(0)) {
                //                        countImgVoxels++;
                //                    }
                //                    if (setForRemoveCuttingIt.GetPixel(p)!=static_cast<InputPixelType>(0)) {
                //                        countRemovingVoxels++;
                //                    }
                //                }

                //Evaluating the removing criteria using the ratio of counting voxels
                //                int ratio = round(( (float)countRemovingVoxels /  (float)countImgVoxels ) * 2.0); //TODO Trocar para o valor de neighborhood certo! - Ver se int arredonda o valor
                //                if(ratio != 0){
                typename InputImageType::SizeType cuttingRadius;
                int statRadiusX = radius[0]+1, statRadiusY = radius[1]+1, statRadiusZ = radius[2]+1; // Increases the calculation area in order to get a more robust statistics
                cuttingRadius[0]=statRadiusX;
                cuttingRadius[1]=statRadiusY;
                cuttingRadius[2]=statRadiusZ;
                //                updateIt.SetRadius(cuttingRadius);
                itk::NeighborhoodIterator<InputImageType, itk::ConstantBoundaryCondition<InputImageType> > inIt(cuttingRadius, auxImage,searchingRegion);
                //                    itk::NeighborhoodIterator<InputImageType, itk::ConstantBoundaryCondition<InputImageType> > inCuttingIt(cuttingRadius, cuttingArea,searchingRegion);

                inIt.SetLocation(imageCuttingIt.GetIndex());
                //                    inCuttingIt.SetLocation(imageCuttingIt.GetIndex());

                /*Calculating the local statistics that will in fact remove the voxel that does not belongs to
                * the brain regions. This estimate is based on the local median values, being here calculated
                * by the increased area of the searching window used for voxel voting process. Hence, each voxel
                * that was previously marked to be investigated as a possible outlier will be or not maintained
                * in the final image by the local median value.
                */
                std::vector<InputPixelType> innerValues;
                unsigned int inN = (statRadiusX*2+1)*(statRadiusY*2+1)*(statRadiusZ*2+1);
                for(unsigned int p = 0; p < inN; p++){
                    innerValues.push_back(inIt.GetPixel(p));
                }
                std::sort(innerValues.begin(), innerValues.end());
                inIt.SetCenterPixel(innerValues[(inN+1)/2]);
                //                }
            }

            ++contourIt;
            ++imageCuttingIt;
            ++setForRemoveCuttingIt;
        }

        // Cleaning the brain mask in order to not propagate this for further iterations
        cuttingArea->FillBuffer(0);

        //Binaring the updated image
        newMask->SetInput(auxImage);
        newMask->SetLowerThreshold(1.0);
        newMask->SetUpperThreshold(itk::NumericTraits<InputPixelType>::max());
        newMask->SetInsideValue(1.0);

        //Choosing the higher binary volume that represents the correct brain volume at iteration i.
        connLabel->SetInput(newMask->GetOutput());
        connLabel->Update();

        relabel->SetInput(connLabel->GetOutput());
        relabel->SetSortByObjectSize(true);
        relabel->Update();

        connectedMask->SetInput(relabel->GetOutput());
        connectedMask->SetLowerThreshold(1);
        connectedMask->SetUpperThreshold(1);
        connectedMask->SetInsideValue(1);

        //Masking the input image with the corrected brain mask
        maskInput->SetInput(updatedImage);
        maskInput->SetMaskImage(connectedMask->GetOutput());
        maskInput->Update();

        /* This is the final iteration step which only update the values for the next iteration.
         * In this case, the image, mask and searching path are updated.
         */
        previousSkeleton->FillBuffer(0);
        itk::ImageRegionIterator<InputImageType>    inIt(maskInput->GetOutput(), maskInput->GetOutput()->GetRequestedRegion());
        itk::ImageRegionIterator<InputImageType>    upIt(updatedImage, updatedImage->GetRequestedRegion());
        itk::ImageRegionIterator<InputImageType>    auxIt(auxImage, auxImage->GetRequestedRegion());
        itk::ImageRegionIterator<LabelImageType>    maskIt(newMask->GetOutput(), newMask->GetOutput()->GetRequestedRegion());
        itk::ImageRegionIterator<LabelImageType>    newMaskIt(newBrainMask, newBrainMask->GetRequestedRegion());
        itk::ImageRegionIterator<InputImageType>    skIt(maskContour->GetOutput(), maskContour->GetOutput()->GetRequestedRegion());
        itk::ImageRegionIterator<InputImageType>    prevSkIt(previousSkeleton, previousSkeleton->GetRequestedRegion());

        upIt.GoToBegin();
        auxIt.GoToBegin();
        inIt.GoToBegin();
        maskIt.GoToBegin();
        newMaskIt.GoToBegin();
        skIt.GoToBegin();
        prevSkIt.GoToBegin();
        while (!upIt.IsAtEnd()) {
            //Copy updated and auxiliary images
            upIt.Set(inIt.Get());
            auxIt.Set(inIt.Get());

            //Copy brain mask
            newMaskIt.Set(maskIt.Get());

            //Copy searching path
            prevSkIt.Set(skIt.Get());

            ++upIt;
            ++auxIt;
            ++inIt;
            ++maskIt;
            ++newMaskIt;
            ++skIt;
            ++prevSkIt;
        }

        //Calculating the brain volume of iteration i+1
        double totalVolumeII;
        itk::ImageRegionIterator<InputImageType>    itII(updatedImage, updatedImage->GetRequestedRegion());
        itII.GoToBegin();
        while (!itII.IsAtEnd()) {
            if (itII.Get()!=static_cast<InputPixelType>(0)) {
                totalVolumeII++;
            }
            ++itII;
        }

        //Stopping criteria
        double c_i = ((totalVolumeI-totalVolumeII)/totalVolumeI);
        std::cout<<"Convergence c = "<<c_i;
        if ( c_i <= convergence) {
            std::cout<<" ( c < "<<convergence<<" )...algorithm converges!"<<std::endl;
            stoppedByConvergence = true;
            break;
        }else{
            std::cout<<" ( c > "<<convergence<<" )...continue iterative process"<<std::endl;
        }
    }

    //After the algorithm being converged, the image and mask (optional) are saved.
    if (updatedMask!="") {
        typedef itk::ImageFileWriter<LabelImageType> WriterType;
        typename WriterType::Pointer labelWriter = WriterType::New();
        labelWriter->SetFileName( updatedMask.c_str() );
        labelWriter->SetInput( newBrainMask );
        labelWriter->SetUseCompression(1);
        labelWriter->Update();
    }
    std::clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    std::cout<<"finished! Total processing time of: "<<elapsed_secs/60.0<<" minutes"<<std::endl;
    std::cout<<"*********************"<<std::endl;

    //Informs what stopping criteria was reached first.
    (stoppedByConvergence)?std::cout<<"Stopped by: Convergence"<<std::endl:std::cout<<"Stopped by: Number of iterations"<<std::endl;

    // Output updated image with less non-brain tissues
    std::cout<<"Output updated image with less non-brain tissues at location: "<<updatedVolume.c_str()<<std::endl;
    if (updatedMask!=""){
        std::cout<<"Output updated mask with less non-brain tissues at location: "<<updatedMask.c_str()<<std::endl;
    }

    //Applying corrected brain mask on the input image
    typename MaskFilterType::Pointer maskingInput = MaskFilterType::New();
    maskingInput->SetInput(inputReader->GetOutput());
    maskingInput->SetMaskImage(newBrainMask);
    maskingInput->Update();

    typedef itk::ImageFileWriter<InputImageType> WriterType;
    typename WriterType::Pointer imageWriter = WriterType::New();
    imageWriter->SetFileName( updatedVolume.c_str() );
    imageWriter->SetInput( maskingInput->GetOutput() );
    imageWriter->SetUseCompression(1);
    imageWriter->Update();

    return EXIT_SUCCESS;
}

} // end of anonymous namespace

int main( int argc, char * argv[] )
{
    PARSE_ARGS;

    itk::ImageIOBase::IOPixelType     pixelType;
    itk::ImageIOBase::IOComponentType componentType;

    try
    {
        itk::GetImageType(inputVolume, pixelType, componentType);

        // This filter handles all types on input, but only produces
        // signed types
        switch( componentType )
        {
        case itk::ImageIOBase::UCHAR:
            return DoIt( argc, argv, static_cast<unsigned char>(0) );
            break;
        case itk::ImageIOBase::CHAR:
            return DoIt( argc, argv, static_cast<signed char>(0) );
            break;
        case itk::ImageIOBase::USHORT:
            return DoIt( argc, argv, static_cast<unsigned short>(0) );
            break;
        case itk::ImageIOBase::SHORT:
            return DoIt( argc, argv, static_cast<short>(0) );
            break;
        case itk::ImageIOBase::UINT:
            return DoIt( argc, argv, static_cast<unsigned int>(0) );
            break;
        case itk::ImageIOBase::INT:
            return DoIt( argc, argv, static_cast<int>(0) );
            break;
        case itk::ImageIOBase::ULONG:
            return DoIt( argc, argv, static_cast<unsigned long>(0) );
            break;
        case itk::ImageIOBase::LONG:
            return DoIt( argc, argv, static_cast<long>(0) );
            break;
        case itk::ImageIOBase::FLOAT:
            return DoIt( argc, argv, static_cast<float>(0) );
            break;
        case itk::ImageIOBase::DOUBLE:
            return DoIt( argc, argv, static_cast<double>(0) );
            break;
        case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
        default:
            std::cerr << "Unknown input image pixel component type: ";
            std::cerr << itk::ImageIOBase::GetComponentTypeAsString( componentType );
            std::cerr << std::endl;
            return EXIT_FAILURE;
            break;
        }
    }

    catch( itk::ExceptionObject & excep )
    {
        std::cerr << argv[0] << ": exception caught !" << std::endl;
        std::cerr << excep << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
