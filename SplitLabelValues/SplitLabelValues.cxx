#include "itkImageFileWriter.h"

#include "itkMaskImageFilter.h"
#include "itkMaskNegatedImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkAddImageFilter.h"

#include "itkPluginUtilities.h"

#include "SplitLabelValuesCLP.h"

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

    typedef unsigned char InputPixelType;
    const unsigned int Dimension = 3;

    typedef itk::Image<InputPixelType,  Dimension> InputImageType;

    typedef itk::ImageFileReader<InputImageType>  ReaderType;

    typename ReaderType::Pointer inputLabelMask = ReaderType::New();
    typename ReaderType::Pointer splitLabelMask = ReaderType::New();

    inputLabelMask->SetFileName( inputLabel.c_str() );
    splitLabelMask->SetFileName( splitLabel.c_str() );

    // Split right hemisphere
    typedef itk::MaskImageFilter<InputImageType, InputImageType>  MaskFilterType;
    typename MaskFilterType::Pointer maskSideA = MaskFilterType::New();
    maskSideA->SetInput(inputLabelMask->GetOutput());
    maskSideA->SetMaskImage(splitLabelMask->GetOutput());
    maskSideA->Update();

    // Change values for the right hemisphere following the rules
    typedef itk::ImageRegionIterator<InputImageType>        RegionIterator;
    RegionIterator      sideAIt(maskSideA->GetOutput(), maskSideA->GetOutput()->GetBufferedRegion());


    if (doKeepSomeValues) {
        bool isSaved;
        sideAIt.GoToBegin();
        while (!sideAIt.IsAtEnd()) {
            isSaved=false;
            for (unsigned int v = 0; v < keepSideA.size(); ++v) {
                if (sideAIt.Get()==keepSideA[v]) {
                    isSaved=true;
                }
            }
            if (!isSaved && sideAIt.Get()!=0) {
                sideAIt.Set(sideAIt.Get()+labelSideA);
            }
            ++sideAIt;
        }
    }else{
        sideAIt.GoToBegin();
        while (!sideAIt.IsAtEnd()) {
            if (sideAIt.Get()!=0) {
                sideAIt.Set(sideAIt.Get()+labelSideA);
            }
            ++sideAIt;
        }
    }

    // Split the left hemispheres
    typedef itk::MaskNegatedImageFilter<InputImageType, InputImageType>   MaskNegatedFilterType;
    typename MaskNegatedFilterType::Pointer maskSideB = MaskNegatedFilterType::New();
    maskSideB->SetInput(inputLabelMask->GetOutput());
    maskSideB->SetMaskImage(splitLabelMask->GetOutput());
    maskSideB->Update();

    // Change values for the right hemisphere following the rules
    typedef itk::ImageRegionIterator<InputImageType>        RegionIterator;
    RegionIterator      sideBIt(maskSideB->GetOutput(), maskSideB->GetOutput()->GetBufferedRegion());

    if (doKeepSomeValues) {
        bool isSaved;
        sideBIt.GoToBegin();
        while (!sideBIt.IsAtEnd()) {
            isSaved=false;
            for (unsigned int v = 0; v < keepSideB.size(); ++v) {
                if (sideBIt.Get()==keepSideB[v]) {
                    isSaved=true;
                }
            }
            if (!isSaved && sideBIt.Get()!=0) {
                sideBIt.Set(sideBIt.Get()+labelSideB);
            }
            ++sideBIt;
        }
    }else{
        sideBIt.GoToBegin();
        while (!sideBIt.IsAtEnd()) {
            if (sideBIt.Get()!=0) {
                sideBIt.Set(sideBIt.Get()+labelSideB);
            }
            ++sideBIt;
        }
    }

    typedef itk::AddImageFilter<InputImageType, InputImageType> AddFilterType;
    typename AddFilterType::Pointer addLabels = AddFilterType::New();
    addLabels->SetInput1(maskSideA->GetOutput());
    addLabels->SetInput2(maskSideB->GetOutput());

    typedef itk::ImageFileWriter<InputImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( outputLabel.c_str() );
    writer->SetInput( addLabels->GetOutput() );
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
