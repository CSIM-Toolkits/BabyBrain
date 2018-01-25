#include "itkImageFileWriter.h"

#include "itkSubtractImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryThinningImageFilter.h"
#include "itkBinaryContourImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkNeighborhoodIterator.h"
#include "itkConstantBoundaryCondition.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"

#include "itkPluginUtilities.h"

#include "BrainVolumeRefinementCLP.h"

#include <ctime>
#include <sstream>
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
    typedef unsigned char LabelPixelType;

    const unsigned int Dimension = 3;

    typedef itk::Image<InputPixelType,  Dimension>    InputImageType;
    typedef itk::Image<LabelPixelType, Dimension>     LabelImageType;

    typedef itk::ImageFileReader<InputImageType>  ImageReaderType;

    typename ImageReaderType::Pointer inputReader = ImageReaderType::New();

    inputReader->SetFileName( inputVolume.c_str() );
    inputReader->Update();

    typedef itk::BinaryThresholdImageFilter<InputImageType, LabelImageType>             BinaryThresholdType;
    typedef itk::BinaryContourImageFilter<LabelImageType, InputImageType>               BinaryContourType;
//    typedef itk::BinaryThinningImageFilter<LabelImageType, LabelImageType>              BinaryThinningType;
    typedef itk::SubtractImageFilter<InputImageType, InputImageType>                    SubtractLabelType;
    typedef itk::GradientMagnitudeImageFilter<InputImageType, InputImageType>           GradientMagnitudeType;
    typedef itk::MedianImageFilter<LabelImageType, LabelImageType>                      MedianFilterType;
    typedef itk::VotingBinaryIterativeHoleFillingImageFilter< LabelImageType >          VotingFillInType;
    typedef itk::ConnectedComponentImageFilter<LabelImageType, LabelImageType>          ConnectedLabelType;
    typedef itk::RelabelComponentImageFilter<LabelImageType, LabelImageType>            RelabelerType;
    typedef itk::BinaryThresholdImageFilter<LabelImageType,LabelImageType>              BinaryLabelMap;
    typedef itk::MaskImageFilter<InputImageType, LabelImageType>                        MaskFilterType;


    // Copying information and pixel values from the input image
    typename InputImageType::Pointer updatedImage = InputImageType::New();
    updatedImage->CopyInformation(inputReader->GetOutput());
    updatedImage->SetRegions(inputReader->GetOutput()->GetRequestedRegion());
    updatedImage->Allocate();
    updatedImage->FillBuffer(0);
    itk::ImageRegionIterator<InputImageType>    inIt(inputReader->GetOutput(), inputReader->GetOutput()->GetRequestedRegion());
    itk::ImageRegionIterator<InputImageType>    upIt(updatedImage, updatedImage->GetRequestedRegion());
    upIt.GoToBegin();
    inIt.GoToBegin();
    while (!upIt.IsAtEnd()) {
        upIt.Set(inIt.Get());
        ++upIt;
        ++inIt;
    }

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

    //Auxiliary skeleton to infer search path for each iteration
    typename InputImageType::Pointer previousSkeleton = InputImageType::New();
    previousSkeleton->CopyInformation(inputReader->GetOutput());
    previousSkeleton->SetRegions(inputReader->GetOutput()->GetRequestedRegion());
    previousSkeleton->Allocate();
    previousSkeleton->FillBuffer(0);

    std::cout<<"*********************"<<std::endl;
    std::cout<<"BrainVolumeRefinement: Initiated"<<std::endl;
    std::cout<<"Iteration: ";
    bool stoppedByConvergence = false;
    for (int n = 0; n < numberOfIterations; ++n) {
        typename BinaryThresholdType::Pointer estimateMask = BinaryThresholdType::New();
        typename BinaryContourType::Pointer maskContour = BinaryContourType::New();
//        typename BinaryThinningType::Pointer skeleton = BinaryThinningType::New();
        typename SubtractLabelType::Pointer actualContour = SubtractLabelType::New();
        typename GradientMagnitudeType::Pointer inputGrad = GradientMagnitudeType::New();
        typename BinaryThresholdType::Pointer newMask = BinaryThresholdType::New();
        typename MedianFilterType::Pointer median = MedianFilterType::New();
        typename VotingFillInType::Pointer fillInHoles = VotingFillInType::New();
        typename ConnectedLabelType::Pointer connLabel = ConnectedLabelType::New();
        typename RelabelerType::Pointer relabel = RelabelerType::New();
        typename BinaryLabelMap::Pointer connectedMask = BinaryLabelMap::New();
        typename MaskFilterType::Pointer maskInput = MaskFilterType::New();

        typename LabelImageType::SizeType radius;
        typename VotingFillInType::InputSizeType vRadius;


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

        //Making mask of the input volume
        std::cout<<(n+1)<<"..."<<std::flush;
        estimateMask->SetInput(updatedImage);
        estimateMask->SetLowerThreshold(1.0);
        estimateMask->SetUpperThreshold(itk::NumericTraits<InputPixelType>::max());
        estimateMask->SetInsideValue( 1 );
        estimateMask->Update();

        // Create volume contour from the input mask
        maskContour->SetInput(estimateMask->GetOutput());
        maskContour->SetForegroundValue( 1 );
        maskContour->Update();

//        skeleton->SetInput(maskContour->GetOutput());
//        skeleton->Update();

//        typedef itk::ImageFileWriter<InputImageType> WriterType3;
//        typename WriterType3::Pointer labelWriter3 = WriterType3::New();
//        std::stringstream out2;
//        out2<<"/home/antonio/Pictures/BVR_test/previousSkeleton"<<n<<"_label.nii.gz";
//        labelWriter3->SetFileName( out2.str() );
//        labelWriter3->SetInput( previousSkeleton );
//        labelWriter3->SetUseCompression(1);
//        labelWriter3->Update();

        //Calculating the actual skeleton for iteration i and its major region
        actualContour->SetInput1(maskContour->GetOutput());
        actualContour->SetInput2(previousSkeleton);
        actualContour->Update();


//        typedef itk::ImageFileWriter<InputImageType> WriterType2;
//        typename WriterType2::Pointer labelWriter2 = WriterType2::New();
//        std::stringstream out;
//        out<<"/home/antonio/Pictures/BVR_test/actualContour"<<n<<"_label.nii.gz";
//        labelWriter2->SetFileName( out.str() );
//        labelWriter2->SetInput( actualContour->GetOutput() );
//        labelWriter2->SetUseCompression(1);
//        labelWriter2->Update();



        //Searching region
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

        // Run through the image gradient (using the countour) in order to find its boundaries
        inputGrad->SetInput(updatedImage);
        inputGrad->Update();

        radius[0] = neighborRadius[0]; // TODO Fazer uma estimativa do tamanho do raio a partir do spacing (1,1,1 pode ser 3,3,3, mas 0.2,0.2,1.0 poderia ser 15,15,3 para manter a mesma proporcao de pixel sendo recrutados... deixar flag bool no XML para isto)
        radius[1] = neighborRadius[1];
        radius[2] = neighborRadius[2];
        itk::ImageRegionIterator<InputImageType> contourIt(actualContour->GetOutput(),searchingRegion);
        itk::NeighborhoodIterator<InputImageType, itk::ConstantBoundaryCondition<InputImageType> > imageIt(radius, updatedImage,searchingRegion);
        itk::NeighborhoodIterator<InputImageType, itk::ConstantBoundaryCondition<InputImageType> > gradIt(radius, inputGrad->GetOutput(),searchingRegion);
        itk::NeighborhoodIterator<InputImageType, itk::ConstantBoundaryCondition<InputImageType> > updateIt(radius, auxImage,searchingRegion);

        // Running over the both image and mask contour in order to find the closest voxel in the brain tissue.
        contourIt.GoToBegin();
        imageIt.GoToBegin();
        gradIt.GoToBegin();
        updateIt.GoToBegin();

        unsigned int N = (neighborRadius[0]*2 + 1)*(neighborRadius[1]*2 + 1)*(neighborRadius[2]*2 + 1);
        InputPixelType wGrad, wIntensity;
        while (!imageIt.IsAtEnd()) {
            if (contourIt.Get()>static_cast<LabelPixelType>(0)) {
                //Calculating statistics into the neighborhood
                wGrad=0.0;
                wIntensity=0.0;
                for (unsigned int p = 0; p < N; p++) {
                    //Weighted average to local gradient
                    wGrad += gradIt.GetPixel(p);
                    //Weighted average to local intensity
                    wIntensity += imageIt.GetPixel(p);
                }
                wGrad/=N;
                wIntensity/=N;

                //Cutting out voxels that does not belongs to the brain based on the mean of the local gradient and mean gray level intensity
                for (unsigned int p = 0; p < N; p++) {
                    if (imageIt.GetPixel(p)!=0) {
                        if (gradIt.GetPixel(p) < (wGrad) && imageIt.GetPixel(p) < (wIntensity)){
                            updateIt.SetPixel(p, 0);
                        }
                    }
                }

            }

            ++contourIt;
            ++imageIt;
            ++updateIt;
            ++gradIt;
        }

        //Binaring the updated image
        newMask->SetInput(auxImage);
        newMask->SetLowerThreshold(1.0);
        newMask->SetUpperThreshold(itk::NumericTraits<InputPixelType>::max());
        newMask->SetInsideValue(1.0);

        if (useMedianFilter) {
            //Median filter in the input image
            median->SetInput(newMask->GetOutput());
            MedianFilterType::RadiusType mRadius;
            mRadius[0] = medianRadius[0];
            mRadius[1] = medianRadius[1];
            mRadius[2] = medianRadius[2];
            median->SetRadius(mRadius);
        }

        //Filling in the holes in the updated brain mask
        //Finding the bigger size of the median filter
        int r = 0;
        for (unsigned int i = 0; i < Dimension; ++i) {
            if (medianRadius[i]>r) {
                r=medianRadius[i];
            }
        }
        vRadius.Fill( r );
        (useMedianFilter)?fillInHoles->SetInput( median->GetOutput() ):fillInHoles->SetInput( newMask->GetOutput() );
        fillInHoles->SetRadius( vRadius );
        fillInHoles->SetMajorityThreshold( majorityThreshold );
        fillInHoles->SetBackgroundValue( 0 );
        fillInHoles->SetForegroundValue( 1 );
        fillInHoles->SetMaximumNumberOfIterations( 100 ); // a default number of iterations to a general case

        //Removing non-connected regions left in the brain mask
        connLabel->SetInput(fillInHoles->GetOutput());
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

        previousSkeleton->FillBuffer(0);
        itk::ImageRegionIterator<InputImageType>    inIt(maskInput->GetOutput(), maskInput->GetOutput()->GetRequestedRegion());
        itk::ImageRegionIterator<InputImageType>    upIt(updatedImage, updatedImage->GetRequestedRegion());
        itk::ImageRegionIterator<InputImageType>    auxIt(auxImage, auxImage->GetRequestedRegion());
        itk::ImageRegionIterator<LabelImageType>    maskIt(connectedMask->GetOutput(), connectedMask->GetOutput()->GetRequestedRegion());
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
            //Copy updated and auxliary images
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
        if ( ((totalVolumeI-totalVolumeII)/totalVolumeI) <= convergence) {
            stoppedByConvergence = true;
            break;
        }
    }
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
    typedef itk::ImageFileWriter<InputImageType> WriterType;
    typename WriterType::Pointer imageWriter = WriterType::New();
    imageWriter->SetFileName( updatedVolume.c_str() );
    imageWriter->SetInput( updatedImage );
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
