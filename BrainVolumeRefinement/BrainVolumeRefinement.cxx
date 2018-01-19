#include "itkImageFileWriter.h"

#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryThinningImageFilter.h"
#include "itkBinaryContourImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodIterator.h"
#include "itkConstantBoundaryCondition.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"

#include "itkPluginUtilities.h"

#include "BrainVolumeRefinementCLP.h"


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
    typedef itk::BinaryContourImageFilter<LabelImageType, LabelImageType>               BinaryContourType;
    typedef itk::BinaryThinningImageFilter<LabelImageType, LabelImageType>              BinaryThinningType;
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

    std::cout<<"*********************"<<std::endl;
    std::cout<<"BrainVolumeRefinement: Initiated"<<std::endl;
    std::cout<<"Iteration: ";
    for (int n = 0; n < numberOfIterations; ++n) {
        typename BinaryThresholdType::Pointer estimateMask = BinaryThresholdType::New();
        typename BinaryContourType::Pointer maskContour = BinaryContourType::New();
        typename BinaryThinningType::Pointer skeleton = BinaryThinningType::New();
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

//        TODO Fazer a subtracao do contorno anterior...assim ele nao passa nas regioes onde ficou estavel
//        TODO Tirar alguns passos do countour iterator (passar para neighborhood iterator) pode ajudar na reducao do tempo...pode colocar uma variavel a mais so para fazer isso
        skeleton->SetInput(maskContour->GetOutput());
        skeleton->Update();

        // Run through the image gradient (using the countour) in order to find its boundaries
        inputGrad->SetInput(updatedImage);
        inputGrad->Update();

        radius[0] = neighborRadius[0];
        radius[1] = neighborRadius[1];
        radius[2] = neighborRadius[2];
        itk::ImageRegionIterator<LabelImageType> contourIt(skeleton->GetOutput(),skeleton->GetOutput()->GetRequestedRegion());
        itk::NeighborhoodIterator<InputImageType, itk::ConstantBoundaryCondition<InputImageType> > imageIt(radius, updatedImage,updatedImage->GetRequestedRegion());
        itk::NeighborhoodIterator<InputImageType, itk::ConstantBoundaryCondition<InputImageType> > gradIt(radius, inputGrad->GetOutput(),inputGrad->GetOutput()->GetRequestedRegion());
        itk::NeighborhoodIterator<InputImageType, itk::ConstantBoundaryCondition<InputImageType> > updateIt(radius, auxImage,auxImage->GetRequestedRegion());

        // Running over the both image and mask contour in order to find the closest voxel in the brain tissue.
        contourIt.GoToBegin();
        imageIt.GoToBegin();
        gradIt.GoToBegin();
        updateIt.GoToBegin();

        unsigned int N = (neighborRadius[0]*2 + 1)*(neighborRadius[1]*2 + 1)*(neighborRadius[2]*2 + 1);
        InputPixelType wGrad, wIntensity;
        while (!imageIt.IsAtEnd()) {
            if (contourIt.Get()!=static_cast<LabelPixelType>(0)) {
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

        //Median filter in the input image
        median->SetInput(newMask->GetOutput());
        MedianFilterType::RadiusType mRadius;
        mRadius[0] = medianRadius[0];
        mRadius[1] = medianRadius[1];
        mRadius[2] = medianRadius[2];
        median->SetRadius(mRadius);

        //Filling in the holes in the updated brain mask
        //Finding the bigger size of the median filter
        int r = 0;
        for (unsigned int i = 0; i < Dimension; ++i) {
            if (medianRadius[i]>r) {
                r=medianRadius[i];
            }
        }
        vRadius.Fill( r );
        fillInHoles->SetInput( median->GetOutput() );
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

        itk::ImageRegionIterator<InputImageType>    inIt(maskInput->GetOutput(), maskInput->GetOutput()->GetRequestedRegion());
        itk::ImageRegionIterator<InputImageType>    upIt(updatedImage, updatedImage->GetRequestedRegion());
        itk::ImageRegionIterator<InputImageType>    auxIt(auxImage, auxImage->GetRequestedRegion());
        upIt.GoToBegin();
        auxIt.GoToBegin();
        inIt.GoToBegin();
        while (!upIt.IsAtEnd()) {
            upIt.Set(inIt.Get());
            auxIt.Set(inIt.Get());
            ++upIt;
            ++auxIt;
            ++inIt;
        }

        if (updatedMask!="" && n == numberOfIterations-1) {

            typedef itk::ImageFileWriter<LabelImageType> WriterType;
            typename WriterType::Pointer labelWriter = WriterType::New();
            labelWriter->SetFileName( updatedMask.c_str() );
            labelWriter->SetInput( connectedMask->GetOutput() );
            labelWriter->SetUseCompression(1);
            labelWriter->Update();
        }

    }
    std::cout<<"finished!"<<std::endl;
    std::cout<<"*********************"<<std::endl;

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
