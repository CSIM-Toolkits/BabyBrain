#include "itkImageFileWriter.h"

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

    typedef itk::BinaryThresholdImageFilter<InputImageType, LabelImageType> BinaryThresholdType;
    typename BinaryThresholdType::Pointer estimateMask = BinaryThresholdType::New();

    //Making mask of the input volume
    std::cout<<"Making mask of the input volume..."<<std::endl;
    estimateMask->SetInput(inputReader->GetOutput());
    estimateMask->SetLowerThreshold(1.0);
    estimateMask->SetUpperThreshold(itk::NumericTraits<InputPixelType>::max());
    estimateMask->SetInsideValue(foregroundValue);

    // Create volume contour from the input mask
    std::cout<<"Creating volume contour from the input mask..."<<std::endl;
    typedef itk::BinaryContourImageFilter<LabelImageType, LabelImageType>   BinaryContourType;
    typename BinaryContourType::Pointer maskContour = BinaryContourType::New();

    maskContour->SetInput(estimateMask->GetOutput());
    maskContour->SetForegroundValue( foregroundValue );
    maskContour->Update();

        typedef itk::BinaryThinningImageFilter<LabelImageType, LabelImageType>  BinaryThinningType;
    typename BinaryThinningType::Pointer skeleton = BinaryThinningType::New();
    skeleton->SetInput(maskContour->GetOutput());
    skeleton->Update();

    // Run through the image gradient (using the countour) in order to find its boundaries
    std::cout<<"Running through the image gradient (using the countour) in order to find its boundaries..."<<std::endl;
    typedef itk::GradientMagnitudeImageFilter<InputImageType, InputImageType> GradientMagnitudeType;
    typename GradientMagnitudeType::Pointer inputGrad = GradientMagnitudeType::New();
    inputGrad->SetInput(inputReader->GetOutput());
    inputGrad->Update();

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

    typename LabelImageType::SizeType radius;
    radius[0] = neighborRadius[0];
    radius[1] = neighborRadius[1];
    radius[2] = neighborRadius[2];
    itk::ImageRegionIterator<LabelImageType> contourIt(skeleton->GetOutput(),skeleton->GetOutput()->GetRequestedRegion());
    itk::NeighborhoodIterator<InputImageType, itk::ConstantBoundaryCondition<InputImageType> > imageIt(radius, inputReader->GetOutput(),inputReader->GetOutput()->GetRequestedRegion());
    itk::NeighborhoodIterator<InputImageType, itk::ConstantBoundaryCondition<InputImageType> > gradIt(radius, inputGrad->GetOutput(),inputGrad->GetOutput()->GetRequestedRegion());
    itk::NeighborhoodIterator<InputImageType, itk::ConstantBoundaryCondition<InputImageType> > updateIt(radius, updatedImage,updatedImage->GetRequestedRegion());

    // Running over the both image and mask contour in order to find the closest voxel in the brain tissue.
    std::cout<<"Running over the both image and mask contour in order to find the closest voxel in the brain tissue..."<<std::endl;
    contourIt.GoToBegin();
    imageIt.GoToBegin();
    gradIt.GoToBegin();
    updateIt.GoToBegin();

    int N = (neighborRadius[0]*2 + 1)*(neighborRadius[1]*2 + 1)*(neighborRadius[2]*2 + 1);
    InputPixelType meanGrad, meanIntensity;
    while (!imageIt.IsAtEnd()) {
        if (contourIt.Get()!=static_cast<LabelPixelType>(0)) {
            //Calculating statistics into the neighborhood
            meanGrad=0.0;
            meanIntensity=0.0;
            for (unsigned int p = 0; p < N; p++) {
                meanGrad += gradIt.GetPixel(p);
                meanIntensity += imageIt.GetPixel(p);
            }
            meanGrad/=N;
            meanIntensity/=N;

            //Cutting out voxels that does not belongs to the brain based on the mean of the local gradient and mean gray level intensity
            for (int p = 0; p < N; p++) {
                if (imageIt.GetPixel(p)!=0) {
                    if (gradIt.GetPixel(p) < (meanGrad) && imageIt.GetPixel(p) < (meanIntensity)){
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
    std::cout<<"Binaring the updated image..."<<std::endl;
    typename BinaryThresholdType::Pointer newMask = BinaryThresholdType::New();
    newMask->SetInput(updatedImage);
    newMask->SetLowerThreshold(1.0);
    newMask->SetUpperThreshold(itk::NumericTraits<InputPixelType>::max());
    newMask->SetInsideValue(1.0);

    //Median filter in the input image
    typedef itk::MedianImageFilter<LabelImageType, LabelImageType>  MedianFilterType;
    typename MedianFilterType::Pointer median = MedianFilterType::New();
    median->SetInput(newMask->GetOutput());
    MedianFilterType::RadiusType mRadius;
    mRadius[0] = medianRadius[0];
    mRadius[1] = medianRadius[1];
    mRadius[2] = medianRadius[2];
    median->SetRadius(mRadius);

    //Filling in the holes in the updated brain mask
    std::cout<<"Filling in the holes in the updated brain mask..."<<std::endl;
    typedef itk::VotingBinaryIterativeHoleFillingImageFilter< LabelImageType > VotingFillInType;
     typename VotingFillInType::InputSizeType vRadius;

    //Finding the bigger size of the median filter
    int r = 0.0;
    for (int i = 0; i < Dimension; ++i) {
        if (medianRadius[i]>r) {
            r=medianRadius[i];
        }
    }
      vRadius.Fill( r );
     typename VotingFillInType::Pointer fillInHoles = VotingFillInType::New();
      fillInHoles->SetInput( median->GetOutput() );
      fillInHoles->SetRadius( vRadius );
      fillInHoles->SetMajorityThreshold( majorityThreshold );
      fillInHoles->SetBackgroundValue( 0 );
      fillInHoles->SetForegroundValue( 1 );
      fillInHoles->SetMaximumNumberOfIterations( numberOfIterations );

    //Masking the input image with the corrected brain mask
    std::cout<<"Masking the input image with the corrected brain mask..."<<std::endl;
    typedef itk::MaskImageFilter<InputImageType, LabelImageType>    MaskFilterType;
    typename MaskFilterType::Pointer maskInput = MaskFilterType::New();
    maskInput->SetInput(inputReader->GetOutput());
    maskInput->SetMaskImage(fillInHoles->GetOutput());
    maskInput->Update();

    // Output updated image with less non-brain tissues
    std::cout<<"Output updated image with less non-brain tissues at location: "<<updatedVolume.c_str()<<std::endl;
    typedef itk::ImageFileWriter<InputImageType> WriterType;
    typename WriterType::Pointer imageWriter = WriterType::New();
    imageWriter->SetFileName( updatedVolume.c_str() );
    imageWriter->SetInput( maskInput->GetOutput() );
    imageWriter->SetUseCompression(1);
    imageWriter->Update();

    if (updatedMask!="") {
        std::cout<<"Output updated mask with less non-brain tissues at location: "<<updatedMask.c_str()<<std::endl;
        typedef itk::ImageFileWriter<LabelImageType> WriterType;
        typename WriterType::Pointer labelWriter = WriterType::New();
        labelWriter->SetFileName( updatedMask.c_str() );
        labelWriter->SetInput( median->GetOutput() );
        labelWriter->SetUseCompression(1);
        labelWriter->Update();
    }

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
