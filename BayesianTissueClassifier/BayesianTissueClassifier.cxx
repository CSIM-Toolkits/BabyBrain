#include "itkImageFileWriter.h"
#include "itkImage.h"

#include "itkPluginUtilities.h"
#include "itkThresholdImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkBayesianClassifierImageFilter.h"
#include "itkBayesianClassifierInitializationImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkImageRegionIterator.h"

#include "BayesianTissueClassifierCLP.h"

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{

template <class T>
int DoIt( int argc, char * argv[], T )
{
    PARSE_ARGS;

    typedef    float                    InputPixelType;
    typedef    unsigned char            OutputPixelType;
    const unsigned int                  Dimension=3;

    typedef itk::Image<InputPixelType,  Dimension> InputImageType;
    typedef itk::Image<OutputPixelType, Dimension> OutputImageType;

    typedef itk::ImageFileReader<InputImageType>    ReaderType;

    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( inputVolume.c_str() );
    reader->Update();

    //Bayesian Segmentation Approach
    typedef itk::BayesianClassifierInitializationImageFilter< InputImageType >         BayesianInitializerType;
    typename BayesianInitializerType::Pointer bayesianInitializer = BayesianInitializerType::New();

    bayesianInitializer->SetInput( reader->GetOutput() );
    bayesianInitializer->SetNumberOfClasses( numberOfTissues );// Background, WM, GM and CSF are defined as 4 tissues
    bayesianInitializer->Update();

    typedef float          PriorType;
    typedef float          PosteriorType;

    typedef itk::VectorImage< InputPixelType, Dimension > VectorInputImageType;
    typedef itk::BayesianClassifierImageFilter< VectorInputImageType,OutputPixelType, PosteriorType,PriorType >   ClassifierFilterType;
    typename ClassifierFilterType::Pointer bayesClassifier = ClassifierFilterType::New();

    bayesClassifier->SetInput( bayesianInitializer->GetOutput() );

    if (imageModality=="T1") {
        if (oneTissue) {
            unsigned char tissueValue = typeTissue;

            //Take only the tissue of interest
            typedef itk::ThresholdImageFilter <OutputImageType>     ThresholdImageFilterType;
            typename  ThresholdImageFilterType::Pointer thresholdFilter  = ThresholdImageFilterType::New();
            thresholdFilter->SetInput(bayesClassifier->GetOutput());
            thresholdFilter->ThresholdOutside(tissueValue, tissueValue);
            thresholdFilter->SetOutsideValue(0);

            //Changing mask value to 1
            typedef itk::DivideImageFilter<OutputImageType, OutputImageType, OutputImageType>   DivideType;
            typename DivideType::Pointer divide = DivideType::New();
            divide->SetInput1(thresholdFilter->GetOutput());
            divide->SetConstant2(tissueValue);
            divide->Update();

            //Removing nonconnected voxels
            typedef unsigned int ConnectedVoxelType;
            typedef itk::Image<ConnectedVoxelType, Dimension>   ConnectedVoxelImageType;
            typedef itk::ConnectedComponentImageFilter<OutputImageType, ConnectedVoxelImageType> ConnectedLabelType;
            typename ConnectedLabelType::Pointer connLabel = ConnectedLabelType::New();
            connLabel->SetInput(divide->GetOutput());
            connLabel->Update();

            typedef itk::RelabelComponentImageFilter<ConnectedVoxelImageType, OutputImageType>      RelabelerType;
            typename RelabelerType::Pointer relabel = RelabelerType::New();
            relabel->SetInput(connLabel->GetOutput());
            relabel->SetSortByObjectSize(true);
            relabel->Update();

            typedef itk::BinaryThresholdImageFilter<OutputImageType,OutputImageType>                BinaryLabelMap;
            typename BinaryLabelMap::Pointer binaryMap = BinaryLabelMap::New();
            binaryMap->SetInput(relabel->GetOutput());
            binaryMap->SetLowerThreshold(1);
            binaryMap->SetUpperThreshold(2);
            binaryMap->SetInsideValue(1);

            typedef itk::ImageFileWriter<OutputImageType> WriterType;
            typename WriterType::Pointer writer = WriterType::New();
            writer->SetFileName( outputLabel.c_str() );
            writer->SetInput( binaryMap->GetOutput() );
            writer->SetUseCompression(1);
            writer->Update();
            return EXIT_SUCCESS;
        }else{
            //        Take all the segmented tissues
            //Removing non connected regions
            typedef unsigned int ConnectedVoxelType;
            typedef itk::Image<ConnectedVoxelType, Dimension>   ConnectedVoxelImageType;
            typedef itk::ConnectedComponentImageFilter<OutputImageType, ConnectedVoxelImageType> ConnectedLabelType;
            typename ConnectedLabelType::Pointer connLabel = ConnectedLabelType::New();
            connLabel->SetInput(bayesClassifier->GetOutput());
            connLabel->Update();

            typedef itk::RelabelComponentImageFilter<ConnectedVoxelImageType, OutputImageType>      RelabelerType;
            typename RelabelerType::Pointer relabel = RelabelerType::New();
            relabel->SetInput(connLabel->GetOutput());
            relabel->SetSortByObjectSize(true);
            relabel->Update();

            typedef itk::ThresholdImageFilter <OutputImageType>     ThresholdImageFilterType;
            typename  ThresholdImageFilterType::Pointer thresholdFilter  = ThresholdImageFilterType::New();
            thresholdFilter->SetInput(relabel->GetOutput());
            thresholdFilter->ThresholdOutside(1, 3);
            thresholdFilter->SetOutsideValue(0);

            typedef itk::MaskImageFilter<OutputImageType, OutputImageType>      MaskTissueType;
            typename MaskTissueType::Pointer maskTissues = MaskTissueType::New();
            maskTissues->SetInput(bayesClassifier->GetOutput());
            maskTissues->SetMaskImage(thresholdFilter->GetOutput());
            maskTissues->SetOutsideValue(0);
            maskTissues->Update();

            //Changing the label values to a list given by the user
            typedef itk::ImageRegionIterator<OutputImageType>        RegionIterator;
            RegionIterator      labelsIt(maskTissues->GetOutput(), maskTissues->GetOutput()->GetRequestedRegion());

            for (int l = 1; l < numberOfTissues; ++l) {
                labelsIt.GoToBegin();
                std::cout<<"Tissue "<<l<<" has value "<<labelsDefinition[l-1]<<std::endl;
                while (!labelsIt.IsAtEnd()) {
                    if (labelsIt.Get()==l) {
                        labelsIt.Set(labelsDefinition[l-1]);
                    }
                    ++labelsIt;
                }
            }

            typedef itk::ImageFileWriter< OutputImageType > WriterType;
            WriterType::Pointer writer = WriterType::New();
            writer->SetInput( maskTissues->GetOutput() );
            writer->SetUseCompression(1);
            writer->SetFileName( outputLabel.c_str() );
            writer->Update();
            return EXIT_SUCCESS;
        }
    }else if (imageModality=="T2") {
        if (oneTissue) {
            unsigned char tissueValue = typeTissue;

            //Take only the tissue of interest
            typedef itk::ThresholdImageFilter <OutputImageType>     ThresholdImageFilterType;
            typename  ThresholdImageFilterType::Pointer thresholdFilter  = ThresholdImageFilterType::New();
            thresholdFilter->SetInput(bayesClassifier->GetOutput());
            thresholdFilter->ThresholdOutside(tissueValue, tissueValue);
            thresholdFilter->SetOutsideValue(0);

            //Changing mask value to 1
            typedef itk::DivideImageFilter<OutputImageType, OutputImageType, OutputImageType>   DivideType;
            typename DivideType::Pointer divide = DivideType::New();
            divide->SetInput1(thresholdFilter->GetOutput());
            divide->SetConstant2(tissueValue);
            divide->Update();

            //Removing nonconnected voxels
            typedef unsigned int ConnectedVoxelType;
            typedef itk::Image<ConnectedVoxelType, Dimension>   ConnectedVoxelImageType;
            typedef itk::ConnectedComponentImageFilter<OutputImageType, ConnectedVoxelImageType> ConnectedLabelType;
            typename ConnectedLabelType::Pointer connLabel = ConnectedLabelType::New();
            connLabel->SetInput(divide->GetOutput());
            connLabel->Update();

            typedef itk::RelabelComponentImageFilter<ConnectedVoxelImageType, OutputImageType>      RelabelerType;
            typename RelabelerType::Pointer relabel = RelabelerType::New();
            relabel->SetInput(connLabel->GetOutput());
            relabel->SetSortByObjectSize(true);
            relabel->Update();

            typedef itk::BinaryThresholdImageFilter<OutputImageType,OutputImageType>                BinaryLabelMap;
            typename BinaryLabelMap::Pointer binaryMap = BinaryLabelMap::New();
            binaryMap->SetInput(relabel->GetOutput());
            binaryMap->SetLowerThreshold(1);
            binaryMap->SetUpperThreshold(2);
            binaryMap->SetInsideValue(1);


            typedef itk::ImageFileWriter<OutputImageType> WriterType;
            typename WriterType::Pointer writer = WriterType::New();
            writer->SetFileName( outputLabel.c_str() );
            writer->SetInput( binaryMap->GetOutput() );
            writer->SetUseCompression(1);
            writer->Update();
            return EXIT_SUCCESS;
        }else{
            //        Take all the segmented tissues
            //Removing non connected regions
            typedef unsigned int ConnectedVoxelType;
            typedef itk::Image<ConnectedVoxelType, Dimension>   ConnectedVoxelImageType;
            typedef itk::ConnectedComponentImageFilter<OutputImageType, ConnectedVoxelImageType> ConnectedLabelType;
            typename ConnectedLabelType::Pointer connLabel = ConnectedLabelType::New();
            connLabel->SetInput(bayesClassifier->GetOutput());
            connLabel->Update();

            typedef itk::RelabelComponentImageFilter<ConnectedVoxelImageType, OutputImageType>      RelabelerType;
            typename RelabelerType::Pointer relabel = RelabelerType::New();
            relabel->SetInput(connLabel->GetOutput());
            relabel->SetSortByObjectSize(true);
            relabel->Update();

            typedef itk::ThresholdImageFilter <OutputImageType>     ThresholdImageFilterType;
            typename  ThresholdImageFilterType::Pointer thresholdFilter  = ThresholdImageFilterType::New();
            thresholdFilter->SetInput(relabel->GetOutput());
            thresholdFilter->ThresholdOutside(1, 3);
            thresholdFilter->SetOutsideValue(0);

            typedef itk::MaskImageFilter<OutputImageType, OutputImageType>      MaskTissueType;
            typename MaskTissueType::Pointer maskTissues = MaskTissueType::New();
            maskTissues->SetInput(bayesClassifier->GetOutput());
            maskTissues->SetMaskImage(thresholdFilter->GetOutput());
            maskTissues->SetOutsideValue(0);
            maskTissues->Update();

            //Changing the label values to a list given by the user
            typedef itk::ImageRegionIterator<OutputImageType>        RegionIterator;
            RegionIterator      labelsIt(maskTissues->GetOutput(), maskTissues->GetOutput()->GetRequestedRegion());

            for (int l = 1; l < numberOfTissues; ++l) {
                labelsIt.GoToBegin();
                std::cout<<"Tissue "<<l<<" has value "<<labelsDefinition[l-1]<<std::endl;
                while (!labelsIt.IsAtEnd()) {
                    if (labelsIt.Get()==l) {
                        labelsIt.Set(labelsDefinition[l-1]);
                    }
                    ++labelsIt;
                }
            }

            typedef itk::ImageFileWriter< OutputImageType > WriterType;
            WriterType::Pointer writer = WriterType::New();
            writer->SetInput( maskTissues->GetOutput() );
            writer->SetUseCompression(1);
            writer->SetFileName( outputLabel.c_str() );
            writer->Update();
            return EXIT_SUCCESS;
        }
    }else if (imageModality=="PD") {
        if (oneTissue) {
            unsigned char tissueValue  = typeTissue;

            //Take only the tissue of interest
            typedef itk::ThresholdImageFilter <OutputImageType>     ThresholdImageFilterType;
            typename  ThresholdImageFilterType::Pointer thresholdFilter  = ThresholdImageFilterType::New();
            thresholdFilter->SetInput(bayesClassifier->GetOutput());
            thresholdFilter->ThresholdOutside(tissueValue, tissueValue);
            thresholdFilter->SetOutsideValue(0);

            //Changing mask value to 1
            typedef itk::DivideImageFilter<OutputImageType, OutputImageType, OutputImageType>   DivideType;
            typename DivideType::Pointer divide = DivideType::New();
            divide->SetInput1(thresholdFilter->GetOutput());
            divide->SetConstant2(tissueValue);
            divide->Update();

            //Removing nonconnected voxels
            typedef unsigned int ConnectedVoxelType;
            typedef itk::Image<ConnectedVoxelType, Dimension>   ConnectedVoxelImageType;
            typedef itk::ConnectedComponentImageFilter<OutputImageType, ConnectedVoxelImageType> ConnectedLabelType;
            typename ConnectedLabelType::Pointer connLabel = ConnectedLabelType::New();
            connLabel->SetInput(divide->GetOutput());
            connLabel->Update();

            typedef itk::RelabelComponentImageFilter<ConnectedVoxelImageType, OutputImageType>      RelabelerType;
            typename RelabelerType::Pointer relabel = RelabelerType::New();
            relabel->SetInput(connLabel->GetOutput());
            relabel->SetSortByObjectSize(true);
            relabel->Update();

            typedef itk::BinaryThresholdImageFilter<OutputImageType,OutputImageType>                BinaryLabelMap;
            typename BinaryLabelMap::Pointer binaryMap = BinaryLabelMap::New();
            binaryMap->SetInput(relabel->GetOutput());
            binaryMap->SetLowerThreshold(1);
            binaryMap->SetUpperThreshold(2);
            binaryMap->SetInsideValue(1);

            typedef itk::ImageFileWriter<OutputImageType> WriterType;
            typename WriterType::Pointer writer = WriterType::New();
            writer->SetFileName( outputLabel.c_str() );
            writer->SetInput( binaryMap->GetOutput() );
            writer->SetUseCompression(1);
            writer->Update();
            return EXIT_SUCCESS;
        }else{
            //        Take all the segmented tissues
            //Removing non connected regions
            typedef unsigned int ConnectedVoxelType;
            typedef itk::Image<ConnectedVoxelType, Dimension>   ConnectedVoxelImageType;
            typedef itk::ConnectedComponentImageFilter<OutputImageType, ConnectedVoxelImageType> ConnectedLabelType;
            typename ConnectedLabelType::Pointer connLabel = ConnectedLabelType::New();
            connLabel->SetInput(bayesClassifier->GetOutput());
            connLabel->Update();

            typedef itk::RelabelComponentImageFilter<ConnectedVoxelImageType, OutputImageType>      RelabelerType;
            typename RelabelerType::Pointer relabel = RelabelerType::New();
            relabel->SetInput(connLabel->GetOutput());
            relabel->SetSortByObjectSize(true);
            relabel->Update();

            typedef itk::ThresholdImageFilter <OutputImageType>     ThresholdImageFilterType;
            typename  ThresholdImageFilterType::Pointer thresholdFilter  = ThresholdImageFilterType::New();
            thresholdFilter->SetInput(relabel->GetOutput());
            thresholdFilter->ThresholdOutside(1, 3);
            thresholdFilter->SetOutsideValue(0);

            typedef itk::MaskImageFilter<OutputImageType, OutputImageType>      MaskTissueType;
            typename MaskTissueType::Pointer maskTissues = MaskTissueType::New();
            maskTissues->SetInput(bayesClassifier->GetOutput());
            maskTissues->SetMaskImage(thresholdFilter->GetOutput());
            maskTissues->SetOutsideValue(0);
            maskTissues->Update();

            //Changing the label values to a list given by the user
            typedef itk::ImageRegionIterator<OutputImageType>        RegionIterator;
            RegionIterator      labelsIt(maskTissues->GetOutput(), maskTissues->GetOutput()->GetRequestedRegion());

            for (int l = 1; l < numberOfTissues; ++l) {
                labelsIt.GoToBegin();
                std::cout<<"Tissue "<<l<<" has value "<<labelsDefinition[l-1]<<std::endl;
                while (!labelsIt.IsAtEnd()) {
                    if (labelsIt.Get()==l) {
                        labelsIt.Set(labelsDefinition[l-1]);
                    }
                    ++labelsIt;
                }
            }

            typedef itk::ImageFileWriter< OutputImageType > WriterType;
            WriterType::Pointer writer = WriterType::New();
            writer->SetInput( maskTissues->GetOutput() );
            writer->SetUseCompression(1);
            writer->SetFileName( outputLabel.c_str() );
            writer->Update();
            return EXIT_SUCCESS;
        }
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
            return DoIt( argc, argv, static_cast<char>(0) );
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
            std::cout << "unknown component type" << std::endl;
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
