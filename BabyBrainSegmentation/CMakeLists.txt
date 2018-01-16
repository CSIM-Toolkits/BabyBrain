#-----------------------------------------------------------------------------
set(MODULE_NAME BabyBrainSegmentation)

#-----------------------------------------------------------------------------
set(MODULE_PYTHON_SCRIPTS
  ${MODULE_NAME}.py
  )

# NEO2012 Brain Atlas
file(GLOB BabyBrain_NEO2012_Atlases_Hemispheres RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "Resources/atlases/NEO2012/*.nii.gz")
file(GLOB BabyBrain_NEO2012_Atlases_Brainmask RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "Resources/atlases/NEO2012/brainmask/*.nii.gz")
file(GLOB BabyBrain_NEO2012_Atlases_Brainstem RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "Resources/atlases/NEO2012/brainstem/*.nii.gz")
file(GLOB BabyBrain_NEO2012_Atlases_Cerebellum RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "Resources/atlases/NEO2012/cerebellum/*.nii.gz")
file(GLOB BabyBrain_NEO2012_Atlases_CSF RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "Resources/atlases/NEO2012/csf/*.nii.gz")
file(GLOB BabyBrain_NEO2012_Atlases_DeepGrayMatter RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "Resources/atlases/NEO2012/dgm/*.nii.gz")
file(GLOB BabyBrain_NEO2012_Atlases_Ventricules RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "Resources/atlases/NEO2012/ventricules/*.nii.gz")
file(GLOB BabyBrain_NEO2012_Templates RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "Resources/atlases/NEO2012/templates/*.nii.gz")

# FET2012 Brain Atlas
file(GLOB BabyBrain_FET2012_Atlases_Hemispheres RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "Resources/atlases/FET2012/*.nii.gz")
file(GLOB BabyBrain_FET2012_Atlases_Brainmask RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "Resources/atlases/FET2012/brainmask/*.nii.gz")
file(GLOB BabyBrain_FET2012_Atlases_Brainstem RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "Resources/atlases/FET2012/brainstem/*.nii.gz")
file(GLOB BabyBrain_FET2012_Atlases_Cerebellum RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "Resources/atlases/FET2012/cerebellum/*.nii.gz")
file(GLOB BabyBrain_FET2012_Atlases_CSF RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "Resources/atlases/FET2012/csf/*.nii.gz")
file(GLOB BabyBrain_FET2012_Atlases_DeepGrayMatter RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "Resources/atlases/FET2012/dgm/*.nii.gz")
file(GLOB BabyBrain_FET2012_Atlases_Ventricules RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "Resources/atlases/FET2012/ventricules/*.nii.gz")
file(GLOB BabyBrain_FET2012_Templates RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "Resources/atlases/FET2012/templates/*.nii.gz")

file(GLOB BabyBrain_ColorTables RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "Resources/*.ctbl")
set(MODULE_PYTHON_RESOURCES
  Resources/Icons/${MODULE_NAME}.png
  ${BabyBrain_NEO2012_Atlases_Hemispheres}
  ${BabyBrain_NEO2012_Atlases_Brainmask}
  ${BabyBrain_NEO2012_Atlases_Brainstem}
  ${BabyBrain_NEO2012_Atlases_Cerebellum}
  ${BabyBrain_NEO2012_Atlases_CSF}
  ${BabyBrain_NEO2012_Atlases_DeepGrayMatter}
  ${BabyBrain_NEO2012_Atlases_Ventricules}
  ${BabyBrain_NEO2012_Templates}

  ${BabyBrain_FET2012_Atlases_Hemispheres}
  ${BabyBrain_FET2012_Atlases_Brainmask}
  ${BabyBrain_FET2012_Atlases_Brainstem}
  ${BabyBrain_FET2012_Atlases_Cerebellum}
  ${BabyBrain_FET2012_Atlases_CSF}
  ${BabyBrain_FET2012_Atlases_DeepGrayMatter}
  ${BabyBrain_FET2012_Atlases_Ventricules}
  ${BabyBrain_FET2012_Templates}

  ${BabyBrain_ColorTables}
  )

#-----------------------------------------------------------------------------
slicerMacroBuildScriptedModule(
  NAME ${MODULE_NAME}
  SCRIPTS ${MODULE_PYTHON_SCRIPTS}
  RESOURCES ${MODULE_PYTHON_RESOURCES}
  WITH_GENERIC_TESTS
  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)

  # Register the unittest subclass in the main script as a ctest.
  # Note that the test will also be available at runtime.
  slicer_add_python_unittest(SCRIPT ${MODULE_NAME}.py)

  # Additional build-time testing
  add_subdirectory(Testing)
endif()
