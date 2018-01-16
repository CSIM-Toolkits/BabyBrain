#include "itkImageFileWriter.h"

#include "itkBinaryThresholdImageFilter.h".h"
#include "itkBinaryContourImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodIterator.h"
#include "itkGradientMagnitudeImageFilter.h"

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
    typedef TPixel OutputPixelType;
    typedef unsigned char LabelPixelType;

    const unsigned int Dimension = 3;

    typedef itk::Image<InputPixelType,  Dimension>    InputImageType;
    typedef itk::Image<OutputPixelType, Dimension>    OutputImageType;
    typedef itk::Image<LabelPixelType, Dimension>     LabelImageType;

    typedef itk::ImageFileReader<InputImageType>  ImageReaderType;
    typedef itk::ImageFileReader<LabelImageType>  LabelReaderType;

    typename ImageReaderType::Pointer inputReader = ImageReaderType::New();
    typename LabelReaderType::Pointer maskReader = LabelReaderType::New();

    inputReader->SetFileName( inputVolume.c_str() );

    typedef itk::BinaryThresholdImageFilter<InputImageType, LabelImageType> BinaryThresholdType;
    typename BinaryThresholdType::Pointer estimateMask = BinaryThresholdType::New();
    if (maskVolume.c_str()=="") {
        //Making mask of the input volume
        estimateMask->SetInput(inputReader->GetOutput());
        estimateMask->SetLowerThreshold(1.0);
        estimateMask->SetUpperThreshold(itk::NumericTraits<InputPixelType>::max());
        estimateMask->SetInsideValue(foregroundValue);
    }else{
        maskReader->SetFileName( maskVolume.c_str() );
    }

    // Create volume contour from the input mask
    typedef itk::BinaryContourImageFilter<LabelImageType, LabelImageType>   BinaryContourType;
    typename BinaryContourType::Pointer maskContour = BinaryContourType::New();

    if (maskVolume.c_str()=="") {
        maskContour->SetInput(estimateMask->GetOutput());
    }else{
        maskContour->SetInput( maskReader->GetOutput() );
    }
    maskContour->SetForegroundValue( foregroundValue );

    // Run through the image gradient (using the countour) in order to find its boundaries
    typedef itk::GradientMagnitudeImageFilter<InputImageType, InputImageType> GradientMagnitudeType;
    typename GradientMagnitudeType::Pointer inputGrad = GradientMagnitudeType::New();
    inputGrad->SetInput(inputReader->GetOutput());

    typename LabelImageType::SizeType radius;
    radius[0] = neighborRadius[0];
    radius[1] = neighborRadius[1];
    radius[2] = neighborRadius[2];
    itk::NeighborhoodIterator<LabelImageType> contourIt(radius, maskContour->GetOutput(),maskContour->GetOutput()->GetRequestedRegion());
    itk::NeighborhoodIterator<InputImageType> imageIt(radius, inputReader->GetOutput(),inputReader->GetOutput()->GetRequestedRegion());
    itk::NeighborhoodIterator<InputImageType> gradIt(radius, inputGrad->GetOutput(),inputGrad->GetOutput()->GetRequestedRegion());


    typename InputImageType::Pointer updatedImage = InputImageType::New();
    updatedImage->CopyInformation(inputReader->GetOutput());
    updatedImage->FillBuffer(0);
    itk::ImageRegionIterator<InputImageType>    pointImageIt(inputReader->GetOutput(), inputReader->GetOutput()->GetRequestedRegion());
    itk::ImageRegionIterator<InputImageType>    pointUpdateIt(updatedImage, updatedImage->GetRequestedRegion());

    pointImageIt.GoToBegin();
    pointUpdateIt.GoToBegin();
    while (!pointImageIt.IsAtEnd()) {
        pointUpdateIt.Set(pointImageIt.Get());
        ++pointImageIt;
        ++pointUpdateIt;
    }

    itk::NeighborhoodIterator<InputImageType> updateIt(radius, updatedImage,updatedImage->GetRequestedRegion());

    // Running over the both image and mask contour in order to find the closest voxel in the brain tissue.
    contourIt.GoToBegin();
    imageIt.GoToBegin();
    updateIt.GoToBegin();
    gradIt.GoToBegin();


    int N = (neighborRadius[0]*2 + 1)*(neighborRadius[1]*2 + 1)*(neighborRadius[2]*2 + 1);
    while (!contourIt.IsAtEnd()) {
        if (contourIt.GetCenterPixel()!=static_cast<LabelPixelType>(0)) {
            //Calculating statistics into the neighborhood
            if (approach=="GradientMagnitude") {
                double mean=0.0;
                for (int p = 0; p < N; ++p) {
                    mean += gradIt.GetPixel(p);
                }
                mean/=N;
                //Cutting out voxels that does not belongs to the brain based on the mean of the local gradient
                for (int p = 0; p < N; ++p) {
                    if (gradIt.GetPixel(p) < (mean)){
                        updateIt.SetPixel(p, 0);
                    }
                }
            }else if (approach=="Mean") {
                double mean=0.0;
                for (int p = 0; p < N; ++p) {
                    mean += imageIt.GetPixel(p);
                }
                mean/=N;
                //Cutting out voxels that does not belongs to the brain based on the mean value
                for (int p = 0; p < N; ++p) {
                    if (imageIt.GetPixel(p) < (mean)){
                        updateIt.SetPixel(p, 0);
                    }
                }
            }

        }
        ++contourIt;
        ++imageIt;
        ++updateIt;
    }

    // Output updated image with less non-brain tissues
    typedef itk::ImageFileWriter<InputImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( updatedVolume.c_str() );
    writer->SetInput( updatedImage );
    writer->SetUseCompression(1);
    writer->Update();

    return EXIT_SUCCESS;
}

//template <typename TPixel>
//void meanNeighbors(typename NeighborhoodIterator neighbors, double mean, double std){
//    N=(neighbors.GetRadius()[0]*2 + 1)*(neighbors.GetRadius()[1]*2 + 1)*(neighbors.GetRadius()[2]*2 + 1);
//    mean=0.0;
//    for (int p = 0; p < N; ++p) {
//        mean += neighbors.GetPixel(p);
//    }
//    mean/=N;

//    std=0.0;
//    for (int p = 0; p < N; ++p) {
//        std += std::pow<double>(mean - neighbors.GetPixel(p),2.0);
//    }
//    std/= std::sqrt(std/(N-1));
//}

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
