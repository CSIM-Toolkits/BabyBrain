#include "itkImageFileWriter.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryContourImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodIterator.h"
#include "itkConstantBoundaryCondition.h"

#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "itkPluginUtilities.h"

#include "LocalLabelingCorrectionCLP.h"

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

    typedef unsigned char LabelPixelType;
    const unsigned int Dimension = 3;

    typedef itk::Image<LabelPixelType,  Dimension> LabelImageType;

    //Read the input label
    typedef itk::ImageFileReader<LabelImageType>  ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( inputLabel.c_str() );
    reader->Update();

    //Creating output label
    LabelImageType::Pointer output = LabelImageType::New();
    output->CopyInformation(reader->GetOutput());
    output->SetRegions(reader->GetOutput()->GetRequestedRegion());
    output->Allocate();
    output->FillBuffer(0);

    //Apply a local decision rule in order to remove wrong segmentation definitions
    typedef itk::BinaryThresholdImageFilter<LabelImageType, LabelImageType>             BinaryThresholdType;
    typename BinaryThresholdType::Pointer thr = BinaryThresholdType::New();
    thr->SetInput(reader->GetOutput());
    thr->SetLowerThreshold( labelToCorrect );
    thr->SetUpperThreshold( labelToCorrect );
    thr->Update();

    //Defining the path where the iterator should run
    typedef itk::BinaryContourImageFilter<LabelImageType, LabelImageType>               BinaryContourType;
    typename BinaryContourType::Pointer maskContour = BinaryContourType::New();
    maskContour->SetInput(thr->GetOutput());
    maskContour->SetForegroundValue( labelToCorrect );
    maskContour->Update();

    typename LabelImageType::SizeType radius;
    radius[0] = neighborRadius[0];
    radius[1] = neighborRadius[1];
    radius[2] = neighborRadius[2];

    typedef itk::NeighborhoodIterator<LabelImageType, itk::ConstantBoundaryCondition<LabelImageType> > NeighborhoodIteratorType;
    typedef itk::ImageRegionIterator<LabelImageType> ImageIteratorType;

    ImageIteratorType           contourIt(maskContour->GetOutput(),maskContour->GetOutput()->GetRequestedRegion());
    ImageIteratorType           outputIt(output, output->GetRequestedRegion());
    NeighborhoodIteratorType    inLabelIt(radius, reader->GetOutput(),reader->GetOutput()->GetRequestedRegion());

    contourIt.GoToBegin();
    inLabelIt.GoToBegin();
    outputIt.GoToBegin();

    unsigned int N = (neighborRadius[0]*2 + 1)*(neighborRadius[1]*2 + 1)*(neighborRadius[2]*2 + 1);
    while (!inLabelIt.IsAtEnd()) {
        if (contourIt.Get()!=static_cast<LabelPixelType>(0)) {
            // Making a local statistics
            double nErrors = 0; //Set the counter for the label errors. This will count how many errors in the local segmentation.
            for (unsigned int p = 0; p < N; p++) {
                if (inLabelIt.GetPixel(p) == static_cast<LabelPixelType>(labelError)){
                    nErrors++;
                }
            }
            nErrors/=N;

            if (nErrors>=tolerance) {
                outputIt.Set(static_cast<LabelPixelType>(labelError));
            }else{
                outputIt.Set(inLabelIt.GetCenterPixel());
            }
        }else{
            //  If the iterator is not on the contour path, the output value is a copy of the input value
            if (inLabelIt.GetCenterPixel()!=static_cast<LabelPixelType>(labelError)) {
                outputIt.Set(inLabelIt.GetCenterPixel());
            }
        }

        ++contourIt;
        ++inLabelIt;
        ++outputIt;
    }


    typedef itk::ImageFileWriter<LabelImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( outputLabel.c_str() );
    writer->SetInput( output );
    writer->SetUseCompression(1);
    writer->Update();

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
        itk::GetImageType(inputLabel, pixelType, componentType);

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
