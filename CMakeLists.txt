cmake_minimum_required(VERSION 3.5)

project(BabyBrain)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "http://slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/BabyBrain")
set(EXTENSION_CATEGORY "Segmentation")
set(EXTENSION_CONTRIBUTORS "Antonio Carlos da S. Senra Filho (University of Sao Paulo)")
set(EXTENSION_DESCRIPTION "This extension ...")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/BabyBrain.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/BabyBrain/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a space separated string, a list or 'NA' if any

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(BabyWholeBrain)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
