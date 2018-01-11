import os
import multiprocessing
import unittest
import vtk, qt, ctk, slicer
import SimpleITK as sitk
import sitkUtils
from slicer.ScriptedLoadableModule import *
import logging

#
# BabyBrainPreparation
#

class BabyBrainPreparation(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "Baby Brain Preparation"
    self.parent.categories = ["Segmentation.Baby Brain"]
    self.parent.dependencies = []
    self.parent.contributors = ["Antonio Carlos da S. Senra Filho (University of Sao Paulo) and Sara"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
This module offers a set of algorithms to biomedical image data preparation, which is focused for the neonate and fetal MRI analysis. The methods used here are
 optimized to structural MRI images, namely T2 and T1, however, any kind of digital 3D images can be processed here (assuming a different set of parameters). 
 The general procedure assumes that the MRI image was already reconstructed in a volume representation and also brain extracted.  
 More details are found in the wikipage: https://www.slicer.org/wiki/Documentation/Nightly/Modules/BabyBrainPreparation
"""
    self.parent.helpText += self.getDefaultModuleDocumentationLink()
    self.parent.acknowledgementText = """
This work was partially funded by CNPq grant ....
""" # replace with organization, grant and thanks.

#
# BabyBrainPreparationWidget
#

class BabyBrainPreparationWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)

    # Instantiate and connect widgets ...

    #
    # Parameters Area
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Input/Output Parameters"
    self.layout.addWidget(parametersCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)

    #
    # input volume selector
    #
    self.inputSelector = slicer.qMRMLNodeComboBox()
    self.inputSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.inputSelector.selectNodeUponCreation = True
    self.inputSelector.addEnabled = False
    self.inputSelector.renameEnabled = True
    self.inputSelector.removeEnabled = True
    self.inputSelector.noneEnabled = False
    self.inputSelector.showHidden = False
    self.inputSelector.showChildNodeTypes = False
    self.inputSelector.setMRMLScene( slicer.mrmlScene )
    self.inputSelector.setToolTip( "Pick the input to the algorithm. This should be an MRI strutural images with a type "
                                   "listed in the Image Modality option." )
    parametersFormLayout.addRow("Input Volume: ", self.inputSelector)

    #
    # output volume selector
    #
    self.outputSelector = slicer.qMRMLNodeComboBox()
    self.outputSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.outputSelector.selectNodeUponCreation = True
    self.outputSelector.addEnabled = True
    self.outputSelector.renameEnabled = True
    self.outputSelector.removeEnabled = True
    self.outputSelector.noneEnabled = False
    self.outputSelector.showHidden = False
    self.outputSelector.showChildNodeTypes = False
    self.outputSelector.setMRMLScene( slicer.mrmlScene )
    self.outputSelector.setToolTip( "Pick the output to the algorithm." )
    parametersFormLayout.addRow("Output Volume: ", self.outputSelector)

    #
    # Enhancement Functions Area
    #
    parametersEnhancementCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersEnhancementCollapsibleButton.text = "Image Enhancement Functions"
    self.layout.addWidget(parametersEnhancementCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersEnhancementLayout = qt.QFormLayout(parametersEnhancementCollapsibleButton)

    #
    # Apply Image Space Resampling
    #
    self.setApplyImageResamplingBooleanWidget = ctk.ctkCheckBox()
    self.setApplyImageResamplingBooleanWidget.setChecked(True)
    self.setApplyImageResamplingBooleanWidget.setToolTip(
      "Apply an image space resampling in the input image in order to reduce the total processing time. This is helpful for "
      "image with large scale space (usually higher than 256^3 mm).")
    parametersEnhancementLayout.addRow("Apply Image Space Resampling", self.setApplyImageResamplingBooleanWidget)

    #
    # Apply AAD filtering
    #
    self.setApplyAADBooleanWidget = ctk.ctkCheckBox()
    self.setApplyAADBooleanWidget.setChecked(True)
    self.setApplyAADBooleanWidget.setToolTip(
      "Apply the AAD filter on the input data. This is recommended because the image noise level may affect the segmentation "
      "performance.")
    parametersEnhancementLayout.addRow("Apply AAD filter", self.setApplyAADBooleanWidget)


    #
    # Apply Bias Correction
    #
    self.setApplyBiasCorrectionBooleanWidget = ctk.ctkCheckBox()
    self.setApplyBiasCorrectionBooleanWidget.setChecked(True)
    self.setApplyBiasCorrectionBooleanWidget.setToolTip(
      "Apply a bias field correction in the input data. This is recommended because the global signal fluctuation provided by "
      "magnetic field inhomogeneity may affect the segmentation performance.")
    parametersEnhancementLayout.addRow("Apply Bias Field Correction", self.setApplyBiasCorrectionBooleanWidget)

    #
    # Apply Global Contrast Enhancement Function
    #
    self.setApplyGlobalEnhancementBooleanWidget = ctk.ctkCheckBox()
    self.setApplyGlobalEnhancementBooleanWidget.setChecked(True)
    self.setApplyGlobalEnhancementBooleanWidget.setToolTip(
      "Apply a full image space contrast enhancement in order to increase tissues delineation. This is recommended because the "
      "original image may present a poor separation between brain tissues.")
    parametersEnhancementLayout.addRow("Apply Global Contrast Enhancement Function", self.setApplyGlobalEnhancementBooleanWidget)


    #
    # Image Space Resampling Parameters Area
    #
    parametersImageResamplingCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersImageResamplingCollapsibleButton.text = "Image Space Resampling Parameters"
    parametersImageResamplingCollapsibleButton.collapsed = True
    self.layout.addWidget(parametersImageResamplingCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersImageResamplingLayout = qt.QFormLayout(parametersImageResamplingCollapsibleButton)

    #
    # Voxel Resampling Size
    #
    self.setVoxelResolutionLineEditWidget = qt.QLineEdit()
    self.setVoxelResolutionLineEditWidget.setText("1,1,1")
    self.setVoxelResolutionLineEditWidget.setToolTip(
      "Voxel resolution to the image resample function. Images with high voxel resolution will strongly increase the processing time.")
    parametersImageResamplingLayout.addRow("Voxel Resampling Resolution", self.setVoxelResolutionLineEditWidget)

    #
    # Interpolation Functions
    #
    self.setInterpolationFunctionComboBoxWidget = ctk.ctkComboBox()
    self.setInterpolationFunctionComboBoxWidget.addItem("bspline")
    self.setInterpolationFunctionComboBoxWidget.addItem("linear")
    self.setInterpolationFunctionComboBoxWidget.addItem("nearestNeighbor")
    self.setInterpolationFunctionComboBoxWidget.addItem("hamming")
    self.setInterpolationFunctionComboBoxWidget.addItem("cosine")
    self.setInterpolationFunctionComboBoxWidget.addItem("welch")
    self.setInterpolationFunctionComboBoxWidget.addItem("lanczos")
    self.setInterpolationFunctionComboBoxWidget.addItem("blackman")
    self.setInterpolationFunctionComboBoxWidget.setToolTip(
      "Interpolation functions.")
    parametersImageResamplingLayout.addRow("Interpolation ", self.setInterpolationFunctionComboBoxWidget)

    #
    # Noise Attenuation Parameters Area
    #
    parametersNoiseAttenuationCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersNoiseAttenuationCollapsibleButton.text = "Noise Attenuation Parameters"
    parametersNoiseAttenuationCollapsibleButton.collapsed = True
    self.layout.addWidget(parametersNoiseAttenuationCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersNoiseAttenuationLayout = qt.QFormLayout(parametersNoiseAttenuationCollapsibleButton)

    #
    # Apply Automatic Conductance Function
    #
    self.setApplyAutoConductanceBooleanWidget = ctk.ctkCheckBox()
    self.setApplyAutoConductanceBooleanWidget.setChecked(True)
    self.setApplyAutoConductanceBooleanWidget.setToolTip(
      "Apply an automatic conductance function in order to adjust the conductance parameters using the input image as source. This is "
      "recommended because the tissues borders are usually difficult to be empirically infered. If this option is checked, the manual "
      "conductance adjustment will be neglected.")
    parametersNoiseAttenuationLayout.addRow("Apply Automatic Conductance Function", self.setApplyAutoConductanceBooleanWidget)

    #
    # Auto Conductance Functions
    #
    self.setConductanceFunctionComboBoxWidget = ctk.ctkComboBox()
    self.setConductanceFunctionComboBoxWidget.addItem("Canny")
    self.setConductanceFunctionComboBoxWidget.addItem("MAD")
    self.setConductanceFunctionComboBoxWidget.addItem("Morphological")
    self.setConductanceFunctionComboBoxWidget.setToolTip(
      "Automatic conductance functions.")
    parametersNoiseAttenuationLayout.addRow("Conductance Function ", self.setConductanceFunctionComboBoxWidget)

    #
    # Apply Conductance Regularization
    #
    self.setApplyConductanceRegularizationBooleanWidget = ctk.ctkCheckBox()
    self.setApplyConductanceRegularizationBooleanWidget.setChecked(True)
    self.setApplyConductanceRegularizationBooleanWidget.setToolTip(
      "Apply an regularization level to the automatic conductance value obtained above. This is recommended only if the automatic conductance "
      "function is overestimating the edge cut level, which may return undesired results (tissue border deterioration). This factor should be "
      "interpreted as a reduction level for less conservative edge preservation. This option will be applyed only if the automatic conductance "
      "choice is used.")
    parametersNoiseAttenuationLayout.addRow("Apply Conductance Regularization", self.setApplyConductanceRegularizationBooleanWidget)

    #
    # Conductance Regularization Value
    #
    self.setConductanceRegularizationIntegerWidget = ctk.ctkSliderWidget()
    self.setConductanceRegularizationIntegerWidget.maximum = 5.0
    self.setConductanceRegularizationIntegerWidget.minimum = 0.1
    self.setConductanceRegularizationIntegerWidget.value = 2.0
    self.setConductanceRegularizationIntegerWidget.singleStep = 0.1
    self.setConductanceRegularizationIntegerWidget.setToolTip("Conductance regularization.")
    parametersNoiseAttenuationLayout.addRow("Conductance Regularization", self.setConductanceRegularizationIntegerWidget)

    #
    # AAD Filtering: Manual Adjustments
    #
    self.groupAADManualParametersBoxButtons = qt.QGroupBox("Manual Adjustments")
    ManualAdjustmentsLayout = qt.QFormLayout()
    self.groupAADManualParametersBoxButtons.setLayout(ManualAdjustmentsLayout)

    self.setFilteringCondutanceWidget = ctk.ctkSliderWidget()
    self.setFilteringCondutanceWidget.maximum = 50
    self.setFilteringCondutanceWidget.minimum = 0
    self.setFilteringCondutanceWidget.value = 10
    self.setFilteringCondutanceWidget.singleStep = 0.1
    self.setFilteringCondutanceWidget.setToolTip("Conductance parameter.")

    self.setFilteringNumberOfIterationWidget = ctk.ctkSliderWidget()
    self.setFilteringNumberOfIterationWidget.maximum = 50
    self.setFilteringNumberOfIterationWidget.minimum = 1
    self.setFilteringNumberOfIterationWidget.value = 25
    self.setFilteringNumberOfIterationWidget.singleStep = 1
    self.setFilteringNumberOfIterationWidget.setToolTip("Number of iterations parameter.")

    self.setFilteringQWidget = ctk.ctkSliderWidget()
    self.setFilteringQWidget.singleStep = 0.01
    self.setFilteringQWidget.minimum = 0.01
    self.setFilteringQWidget.maximum = 1.99
    self.setFilteringQWidget.value = 1.4
    self.setFilteringQWidget.setToolTip("Q value parameter.")

    self.setFilteringTimeStepWidget = ctk.ctkSliderWidget()
    self.setFilteringTimeStepWidget.singleStep = 0.0001
    self.setFilteringTimeStepWidget.minimum = 0.0001
    self.setFilteringTimeStepWidget.maximum = 0.0625
    self.setFilteringTimeStepWidget.value = 0.04
    self.setFilteringTimeStepWidget.setToolTip("Time step parameter.")

    ManualAdjustmentsLayout.addRow("Conductance ", self.setFilteringCondutanceWidget)
    ManualAdjustmentsLayout.addRow("Number Of Iterations ", self.setFilteringNumberOfIterationWidget)
    ManualAdjustmentsLayout.addRow("Q Value ", self.setFilteringQWidget)
    ManualAdjustmentsLayout.addRow("Time Step ", self.setFilteringTimeStepWidget)

    parametersNoiseAttenuationLayout.addRow(self.groupAADManualParametersBoxButtons)

    #
    # Bias Field Correction Parameters Area
    #
    parametersBiasCorrectionCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersBiasCorrectionCollapsibleButton.text = "Bias Field Correction Parameters"
    parametersBiasCorrectionCollapsibleButton.collapsed = True
    self.layout.addWidget(parametersBiasCorrectionCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersBiasCorrectionLayout = qt.QFormLayout(parametersBiasCorrectionCollapsibleButton)

    #
    # input mask selector
    #
    self.inputMaskSelector = slicer.qMRMLNodeComboBox()
    self.inputMaskSelector.nodeTypes = ["vtkMRMLLabelMapVolumeNode"]
    self.inputMaskSelector.selectNodeUponCreation = True
    self.inputMaskSelector.addEnabled = False
    self.inputMaskSelector.removeEnabled = True
    self.inputMaskSelector.noneEnabled = True
    self.inputMaskSelector.showHidden = False
    self.inputMaskSelector.showChildNodeTypes = False
    self.inputMaskSelector.setMRMLScene( slicer.mrmlScene )
    self.inputMaskSelector.setToolTip( "Pick the input brain mask. This image mask should informs the region of interest "
                                       "where the bias field correction will be applied. If not provided, an Otsu segmentation"
                                       " is performed." )
    parametersBiasCorrectionLayout.addRow("Image Mask: ", self.inputMaskSelector)

    #
    # Spline Grid - First Level
    #
    self.setSplineGridFirstLevelLineEditWidget = qt.QLineEdit()
    self.setSplineGridFirstLevelLineEditWidget.setText("1,1,1")
    self.setSplineGridFirstLevelLineEditWidget.setToolTip(
      "Spline grid for bias field sampling. This process is separated in two levels using one before and other after the noise attenuation process. "
      "This parameter is supposed to be conjugated with the Shrink Factor, i.e., the higher the spline grid the lower should be the Shrink Factor.")
    parametersBiasCorrectionLayout.addRow("Spline Grid (First Level) ", self.setSplineGridFirstLevelLineEditWidget)

    #
    # Spline Grid - Second Level
    #
    self.setSplineGridSecondLevelLineEditWidget = qt.QLineEdit()
    self.setSplineGridSecondLevelLineEditWidget.setText("4,4,4")
    self.setSplineGridSecondLevelLineEditWidget.setToolTip(
      "The same as in Spline Grid (First Level).")
    parametersBiasCorrectionLayout.addRow("Spline Grid (Second Level) ", self.setSplineGridSecondLevelLineEditWidget)

    #
    # Shrink Factor - First Level
    #
    self.setShrinkFactorFirstLevelWidget = qt.QSpinBox()
    self.setShrinkFactorFirstLevelWidget.setMinimum(1)
    self.setShrinkFactorFirstLevelWidget.setMaximum(30)
    self.setShrinkFactorFirstLevelWidget.setValue(4)
    self.setShrinkFactorFirstLevelWidget.setToolTip(
      "Defines how much the image should be upsampled before estimating the inhomogeneity field. Increase if you want to reduce the execution"
      " time. 1 corresponds to the original resolution. Larger values will significantly reduce the computation time, however may not resolve"
      "images with strong field inhomogeneity.")
    parametersBiasCorrectionLayout.addRow("Shrink Factor (First Level) ", self.setShrinkFactorFirstLevelWidget)

    #
    # Shrink Factor - Second Level
    #
    self.setShrinkFactorSecondLevelWidget = qt.QSpinBox()
    self.setShrinkFactorSecondLevelWidget.setMinimum(1)
    self.setShrinkFactorSecondLevelWidget.setMaximum(30)
    self.setShrinkFactorSecondLevelWidget.setValue(3)
    self.setShrinkFactorSecondLevelWidget.setToolTip(
      "The same as in Shrink Factor (First Level).")
    parametersBiasCorrectionLayout.addRow("Shrink Factor (Second Level) ", self.setShrinkFactorSecondLevelWidget)

    #
    # Global Contrast Enhancement Parameters Area
    #
    parametersGlobalContrastEnhancementCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersGlobalContrastEnhancementCollapsibleButton.text = "Global Contrast Enhancement Parameters"
    parametersGlobalContrastEnhancementCollapsibleButton.collapsed = True
    self.layout.addWidget(parametersGlobalContrastEnhancementCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersGlobalContrastEnhancementLayout = qt.QFormLayout(parametersGlobalContrastEnhancementCollapsibleButton)

    #
    # Contrast Modulation Functions
    #
    self.setContrastModulationFunctionComboBoxWidget = ctk.ctkComboBox()
    self.setContrastModulationFunctionComboBoxWidget.addItem("Logistic")
    self.setContrastModulationFunctionComboBoxWidget.addItem("CLAHE") #TODO Implementar CLAHE e CDF image contrast enhancement functions para ITK!!!
    self.setContrastModulationFunctionComboBoxWidget.addItem("CDF")
    self.setContrastModulationFunctionComboBoxWidget.setToolTip(
      "Contrast modulation functions used to enhance tissue signal contrast in all image scale space. These methods are based on image histogram.")
    parametersGlobalContrastEnhancementLayout.addRow("Contrast Modulation Function ", self.setContrastModulationFunctionComboBoxWidget)

    #
    # Contrast Modulation: Logistic Parameters
    #
    self.groupGlobalContrastParametersBoxButtons = qt.QGroupBox("Logistic Parameters")
    LogisticParametersLayout = qt.QFormLayout()
    self.groupGlobalContrastParametersBoxButtons.setLayout(LogisticParametersLayout)

    self.setHigherCutWidget = ctk.ctkSliderWidget() # TODO Ver como fazer um doubleSpingBox parecido com o windowsLevel...
    self.setHigherCutWidget.maximum = 0.99
    self.setHigherCutWidget.minimum = 0.01
    self.setHigherCutWidget.value = 0.96
    self.setHigherCutWidget.singleStep = 0.01
    self.setHigherCutWidget.setToolTip("Higher cut parameter.")

    self.setLowerCutWidget = ctk.ctkSliderWidget()
    self.setLowerCutWidget.maximum = 0.99
    self.setLowerCutWidget.minimum = 0.01
    self.setLowerCutWidget.value = 0.04
    self.setLowerCutWidget.singleStep = 0.01
    self.setLowerCutWidget.setToolTip("Lower cut parameter.")

    self.setMaximumOutputWeightingWidget = qt.QSpinBox()
    self.setMaximumOutputWeightingWidget.setMinimum(1)
    # self.setMaximumOutputWeightingWidget.setMaximum(30)
    self.setMaximumOutputWeightingWidget.setValue(2)
    self.setMaximumOutputWeightingWidget.setToolTip(
      "This informs the maximum scaling factor to the contrast enhancement method. The higher it is, the stronger will be contrast modulation.")

    self.setMinimumOutputWeightingWidget = qt.QSpinBox()
    self.setMinimumOutputWeightingWidget.setMinimum(1)
    # self.setMinimumOutputWeightingWidget.setMaximum(30)
    self.setMinimumOutputWeightingWidget.setValue(1)
    self.setMinimumOutputWeightingWidget.setToolTip(
      "This informs the minimum scaling factor to the contrast enhancement method. The lower it is, the stronger will be contrast modulation.")

    self.setFlipWeightingFunctionBooleanWidget = ctk.ctkCheckBox()
    self.setFlipWeightingFunctionBooleanWidget.setChecked(False)
    self.setFlipWeightingFunctionBooleanWidget.setToolTip(
      "Check this if you want to flip the contrast enhancement calculation in the input image. This will invert the weighting function.")

    LogisticParametersLayout.addRow("Higher Cut ", self.setHigherCutWidget)
    LogisticParametersLayout.addRow("Lower Cut ", self.setLowerCutWidget)
    LogisticParametersLayout.addRow("Maximum Weighting ", self.setMaximumOutputWeightingWidget)
    LogisticParametersLayout.addRow("Minimum Weighting ", self.setMinimumOutputWeightingWidget)
    LogisticParametersLayout.addRow("Flip Weighting Function ", self.setFlipWeightingFunctionBooleanWidget)

    parametersGlobalContrastEnhancementLayout.addRow(self.groupGlobalContrastParametersBoxButtons)

    #
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Apply")
    self.applyButton.toolTip = "Run the algorithm."
    self.applyButton.enabled = False
    parametersFormLayout.addRow(self.applyButton)

    # connections
    self.applyButton.connect('clicked(bool)', self.onApplyButton)
    self.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.outputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSelect()

  def cleanup(self):
    pass

  def onSelect(self):
    self.applyButton.enabled = self.inputSelector.currentNode() and self.outputSelector.currentNode()


  def onApplyButton(self):
    logic = BabyBrainPreparationLogic()
    # imageModality = self.setImageModalityBooleanWidget.currentText
    useResampling = self.setApplyImageResamplingBooleanWidget.isChecked()
    useAAD = self.setApplyAADBooleanWidget.isChecked()
    useN4ITK = self.setApplyBiasCorrectionBooleanWidget.isChecked()
    useGlobalContrastEnhancement = self.setApplyGlobalEnhancementBooleanWidget.isChecked()
    voxelResampling = self.setVoxelResolutionLineEditWidget.text
    interpolation = self.setInterpolationFunctionComboBoxWidget.currentText
    applyAutoConductance = self.setApplyAutoConductanceBooleanWidget.isChecked()
    conductanceFunction = self.setConductanceFunctionComboBoxWidget.currentText
    applyConductanceRegularization = self.setApplyConductanceRegularizationBooleanWidget.isChecked()
    conductanceRegularization = self.setConductanceRegularizationIntegerWidget.value
    conductance = self.setFilteringCondutanceWidget.value
    numInt = self.setFilteringNumberOfIterationWidget.value
    qValue = self.setFilteringQWidget.value
    timeStep = self.setFilteringTimeStepWidget.value
    splineGrid_1 = self.setSplineGridFirstLevelLineEditWidget.text
    splineGrid_2 = self.setSplineGridSecondLevelLineEditWidget.text
    shrinkFactor_1 = self.setShrinkFactorFirstLevelWidget.value
    shrinkFactor_2 = self.setShrinkFactorSecondLevelWidget.value
    globalContrastEnhancement = self.setContrastModulationFunctionComboBoxWidget.currentText
    higherCut = self.setHigherCutWidget.value
    lowerCut = self.setLowerCutWidget.value
    maxWeight = self.setMaximumOutputWeightingWidget.value
    minWeight = self.setMinimumOutputWeightingWidget.value
    useFlip = self.setFlipWeightingFunctionBooleanWidget.isChecked()
    logic.run(self.inputSelector.currentNode()
              , self.outputSelector.currentNode()
              , useResampling
              , useAAD
              , useN4ITK
              , useGlobalContrastEnhancement
              , voxelResampling
              , interpolation
              , applyAutoConductance
              , conductanceFunction
              , applyConductanceRegularization
              , conductanceRegularization
              , conductance
              , numInt
              , qValue
              , timeStep
              , self.inputMaskSelector.currentNode()
              , splineGrid_1
              , splineGrid_2
              , shrinkFactor_1
              , shrinkFactor_2
              , globalContrastEnhancement
              , higherCut
              , lowerCut
              , maxWeight
              , minWeight
              , useFlip)

#
# BabyBrainPreparationLogic
#

class BabyBrainPreparationLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """
  def hasImageData(self,volumeNode):
    """This is an example logic method that
    returns true if the passed in volume
    node has valid image data
    """
    if not volumeNode:
      logging.debug('hasImageData failed: no volume node')
      return False
    if volumeNode.GetImageData() is None:
      logging.debug('hasImageData failed: no image data in volume node')
      return False
    return True

  def isValidInputOutputData(self, inputVolumeNode, outputVolumeNode):
    """Validates if the output is not the same as input
    """
    if not inputVolumeNode:
      logging.debug('isValidInputOutputData failed: no input volume node defined')
      return False
    if not outputVolumeNode:
      logging.debug('isValidInputOutputData failed: no output volume node defined')
      return False
    if inputVolumeNode.GetID()==outputVolumeNode.GetID():
      logging.debug('isValidInputOutputData failed: input and output volume is the same. Create a new volume for output to avoid this error.')
      return False
    return True


  def run(self, inputVolume
          , outputVolume
          , useResampling
          , useAAD
          , useN4ITK
          , useGlobalContrastEnhancement
          , voxelResampling
          , interpolation
          , applyAutoConductance
          , conductanceFunction
          , applyConductanceRegularization
          , conductanceRegularization
          , conductance
          , numInt
          , qValue
          , timeStep
          , maskVolume
          , splineGrid_1
          , splineGrid_2
          , shrinkFactor_1
          , shrinkFactor_2
          , globalContrastEnhancement
          , higherCut
          , lowerCut
          , maxWeight
          , minWeight
          , useFlip):
    """
    Run the actual algorithm
    """

    if not self.isValidInputOutputData(inputVolume, outputVolume):
      slicer.util.errorDisplay('Input volume is the same as output volume. Choose a different output volume.')
      return False

    logging.info('Processing started')

    # Copying the input node into the output node in order to apply the changes in the image and not changing the input data
    volumesLogic = slicer.modules.volumes.logic()
    volumesLogic.CreateScalarVolumeFromVolume(slicer.mrmlScene, outputVolume, inputVolume)

    # TODO Tentar novamente criar uma barra de progresso : slicer.util.createProgressDialog
    # progressbar = slicer.util.createProgressDialog(parent=slicer.util.mainWindow(),autoClose=False)

    # TODO Verificar se fazer a separacao brainstem+cerebellum do resto pode ser melhor para poder filtrar com o filtro AAD. A imagem enh.nrrd parece q tem distorce as bordas do cerebello...dai fica ruim a segmentao da fossa posterior

    if useResampling:
        #
        # Step 1 - Space Resampling
        #
        self.imageResamplingResolution(inputVolume
                                       , outputVolume
                                       , voxelResampling
                                       , interpolation)
        slicer.util.showStatusMessage("Image resampling space is finished...")
    # slicer.app.processEvents()
    # progressbar.value = 20
    # progressbar.labelText = "Space Resampling..."

    if useN4ITK:
        #
        # Step 2 - Bias field correction - level 1
        #
        self.biasFielCorretion(outputVolume
                               , outputVolume
                               , maskVolume
                               , splineGrid_1
                               , shrinkFactor_1
                               , useResampling
                               , voxelResampling)
        slicer.util.showStatusMessage("Bias field correction (first level) is finished...")

    if useAAD:
        #
        # Step 3 - Noise attenuation
        #
        self.anisotropicAnomalousDiffusionFilter(outputVolume
                                                 , outputVolume
                                                 , numInt
                                                 , qValue
                                                 , timeStep
                                                 , conductance
                                                 , applyAutoConductance
                                                 , conductanceFunction
                                                 , applyConductanceRegularization
                                                 , conductanceRegularization)
        slicer.util.showStatusMessage("Image noise attenuation is finished...")
        # progressbar.value = 40
        # progressbar.labelText = "Image noise attenuation..."

    if useN4ITK:
        #
        # Step  - Bias field correction - level 2
        #
        self.biasFielCorretion(outputVolume
                               , outputVolume
                               , maskVolume
                               , splineGrid_2
                               , shrinkFactor_2
                               , useResampling
                               , voxelResampling)
        slicer.util.showStatusMessage("Bias field correction (second level) is finished...")

    if useGlobalContrastEnhancement:
        #
        # Step  - Global image contrast enhancement
        #
        self.imageGlobalConstrastEnhancement(outputVolume
                                       , outputVolume
                                       , globalContrastEnhancement
                                       , lowerCut
                                       , higherCut
                                       , maxWeight
                                       , minWeight
                                       , useFlip)
        slicer.util.showStatusMessage("Global contrast enhancement is finished...")

    # progressbar.value = 100
    # slicer.app.processEvents()
    # progressbar.labelText = "Baby Brain - Preparation is finished!"

    logging.info('Processing completed')
    slicer.util.showStatusMessage("Baby Brain Preparation is finished!")

    return True

  #
  # Image Resampling Resolution
  #
  def imageResamplingResolution(self, inputNode
                        , outputNode
                        , resolution
                        , interpolation):
    params = {}
    params["InputVolume"] = inputNode.GetID()
    params["OutputVolume"] = outputNode.GetID()
    params["outputPixelSpacing"] = resolution
    params["interpolationType"] = interpolation

    slicer.cli.run(slicer.modules.resamplescalarvolume, None, params, wait_for_completion=True)

  #
  # N4ITK Bias Field Correction
  #
  def biasFielCorretion(self, inputNode
                        , outputNode
                        , maskNode
                        , splineGrid
                        , shrinkFactor
                        , useResampling
                        , voxelResampling):
    params = {}
    params["inputImageName"] = inputNode.GetID()
    params["outputImageName"] = outputNode.GetID()
    if maskNode is not None:
        if useResampling:
            #Resampling the input mask
            self.imageResamplingResolution(maskNode,maskNode,voxelResampling,"nearestNeighbor")
            params["maskImageName"] = maskNode.GetID()
        else:
            params["maskImageName"] = maskNode.GetID()
    params["initialMeshResolution"] = splineGrid
    params["shrinkFactor"] = shrinkFactor

    slicer.cli.run(slicer.modules.n4itkbiasfieldcorrection, None, params, wait_for_completion=True)

  #
  # Anisotropic Anomalous Diffusion Filter
  #
  def anisotropicAnomalousDiffusionFilter(self, inputNode
                                          , outputNode
                                          , numberOfIterations
                                          , qValue
                                          , timeStep
                                          , conductance
                                          , useAutoConductance
                                          , conductanceFunction
                                          , useConductanceRegularization
                                          , conductanceRegularization):
    params = {}
    params["inputVolume"] = inputNode.GetID()
    params["outputVolume"] = outputNode.GetID()
    params["useAutoConductance"] = useAutoConductance
    params["conductanceFunction"] = conductanceFunction
    params["useConductanceRegularization"] = useConductanceRegularization
    params["kappaRegFactor"] = conductanceRegularization
    params["optFunction"] = "Canny"
    params["iterations"] = numberOfIterations
    params["conductance"] = conductance
    params["q"] = qValue
    params["timeStep"] = timeStep

    slicer.cli.run(slicer.modules.aadimagefilter, None, params, wait_for_completion=True)

  #
  # Global Contrast Enhancement
  #
  def imageGlobalConstrastEnhancement(self, inputNode
                                      , outputNode
                                      , algorithm
                                      , lowerCut
                                      , higherCut
                                      , maximumScaling
                                      , minimumScaling
                                      , flipFunction):
    params = {}
    params["inputVolume"] = inputNode.GetID()
    params["outputVolume"] = outputNode.GetID()
    params["algorithm"] = algorithm
    params["lowerCut"] = lowerCut
    params["higherCut"] = higherCut
    params["maximumScaling"] = maximumScaling
    params["minimumScaling"] = minimumScaling
    params["flipFunction"] = flipFunction

    slicer.cli.run(slicer.modules.globalcontrastenhancer, None, params, wait_for_completion=True)



class BabyBrainPreparationTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_BabyBrainPreparation1()

  def test_BabyBrainPreparation1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")
    #
    # first, get some data
    #
    import urllib
    downloads = (
        ('http://slicer.kitware.com/midas3/download?items=5767', 'FA.nrrd', slicer.util.loadVolume),
        )

    for url,name,loader in downloads:
      filePath = slicer.app.temporaryPath + '/' + name
      if not os.path.exists(filePath) or os.stat(filePath).st_size == 0:
        logging.info('Requesting download %s from %s...\n' % (name, url))
        urllib.urlretrieve(url, filePath)
      if loader:
        logging.info('Loading %s...' % (name,))
        loader(filePath)
    self.delayDisplay('Finished with download and loading')

    volumeNode = slicer.util.getNode(pattern="FA")
    logic = BabyBrainPreparationLogic()
    self.assertIsNotNone( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')
