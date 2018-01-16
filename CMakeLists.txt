cmake_minimum_required(VERSION 3.5)

project(BabyBrain)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "http://slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/BabyBrain")
set(EXTENSION_CATEGORY "Segmentation")
set(EXTENSION_CONTRIBUTORS "Antonio Carlos da S. Senra Filho (University of Sao Paulo)")
set(EXTENSION_DESCRIPTION "This extension offers a set of algorithms to biomedical image data preparation and segmentation, which are focused for the neonate and fetal MRI analysis. The methods used here are
 optimized to structural MRI images, namely T2 and T1, however, any kind of digital 3D images can be processed here (assuming a different set of parameters). ...ADICIONAR UMA EXPLICACAO DOS METODOS...
 The general procedure assumes that the MRI image was already reconstructed in a volume representation and also brain extracted.")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/BabyBrain.png")
set(EXTENSION_SCREENSHOTURLS "http://.../Extensions/BabyBrain/Screenshots/1.png")
set(EXTENSION_DEPENDS "AnomalousFiltersExtension") # Specified as a space separated string, a list or 'NA' if any

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(BabyBrainSegmentation)
add_subdirectory(BabyBrainPreparation)
add_subdirectory(GlobalContrastEnhancer)
add_subdirectory(BayesianTissueClassifier)
add_subdirectory(SplitLabelValues)
add_subdirectory(BrainVolumeRefinement)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
