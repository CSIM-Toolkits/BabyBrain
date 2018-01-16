import os
import multiprocessing
import platform
import sys
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging

import SimpleITK as sitk
import sitkUtils

from os.path import expanduser

#
# BabyBrainSegmentation
#

class BabyBrainSegmentation(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "Baby Brain Segmentation"
    self.parent.categories = ["Segmentation.Baby Brain"]
    self.parent.dependencies = []
    self.parent.contributors = ["Antonio Carlos da S. Senra Filho (University of Sao Paulo) and Sara"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
This module offers a brain tissue segmentation pipeline to fetal, neonatal and pediatric MRI images. At moment, the structural MRI images are supported
, namely T1, T2 and PD. The general segmentation sequence is based on a naive Bayes classifier coupled to a local signal smoothing and label
propagation step from an brain atlas. More details are found at the wikipage: https://www.slicer.org/wiki/Documentation/Nightly/Modules/BabyBrainSegmentation 
"""
    self.parent.helpText += self.getDefaultModuleDocumentationLink()
    self.parent.acknowledgementText = """
This work was partially funded by CNPq grant ....
""" # replace with organization, grant and thanks.

#
# BabyBrainSegmentationWidget
#

class BabyBrainSegmentationWidget(ScriptedLoadableModuleWidget):
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
    self.inputSelector.removeEnabled = True
    self.inputSelector.renameEnabled = True
    self.inputSelector.noneEnabled = False
    self.inputSelector.showHidden = False
    self.inputSelector.showChildNodeTypes = False
    self.inputSelector.setMRMLScene( slicer.mrmlScene )
    self.inputSelector.setToolTip( "Pick the input to the algorithm. It is recommended using a preprocessed image here (e.g. Baby Brain Preparation module output)" )
    parametersFormLayout.addRow("Input Volume: ", self.inputSelector)

    #
    # Image Modality
    #
    self.setImageModalityBooleanWidget = ctk.ctkComboBox()
    self.setImageModalityBooleanWidget.addItem("T2")
    self.setImageModalityBooleanWidget.addItem("T1")
    self.setImageModalityBooleanWidget.setToolTip(
      "MRI strutural image inserted as a input volume.")
    parametersFormLayout.addRow("Image Modality ", self.setImageModalityBooleanWidget)

    #
    # Brain tissues selector
    #
    self.brainTissuesSelector = slicer.qMRMLNodeComboBox()
    self.brainTissuesSelector.nodeTypes = ["vtkMRMLLabelMapVolumeNode"]
    self.brainTissuesSelector.selectNodeUponCreation = True
    self.brainTissuesSelector.addEnabled = True
    self.brainTissuesSelector.renameEnabled = True
    self.brainTissuesSelector.removeEnabled = True
    self.brainTissuesSelector.noneEnabled = False
    self.brainTissuesSelector.showHidden = False
    self.brainTissuesSelector.showChildNodeTypes = False
    self.brainTissuesSelector.setMRMLScene( slicer.mrmlScene )
    self.brainTissuesSelector.setToolTip( "Pick the output brain tissues label." )
    parametersFormLayout.addRow("Brain Tissues: ", self.brainTissuesSelector)

    #
    # Estimate Basal Ganglia Structures?
    #
    self.setUseBasalGangliaEstimatorBooleanWidget = ctk.ctkCheckBox()
    self.setUseBasalGangliaEstimatorBooleanWidget.setChecked(False)
    self.setUseBasalGangliaEstimatorBooleanWidget.setToolTip(
      "Check this if you want to output the global basal ganglia structures in the final segmentation. This step is in "
      "experimental stage which is only based on a label propagation from the chosen brain atlas.")
    parametersFormLayout.addRow("Estimate Basal Ganglia Structures",
                                            self.setUseBasalGangliaEstimatorBooleanWidget)


    # #
    # # Debug Mode
    # #
    # self.setUseDebugModeBooleanWidget = ctk.ctkCheckBox()
    # self.setUseDebugModeBooleanWidget.setChecked(False)
    # self.setUseDebugModeBooleanWidget.setToolTip(
    #   "Check this if you want to keep the intermediate files generated through the segmentation process. This could be"
    #   "useful if you want to check each brain structure labels in a separated manner.")
    # parametersFormLayout.addRow("Debug Mode",
    #                                         self.setUseDebugModeBooleanWidget)

    #
    # Atlas Propagation Parameters Area
    #
    parametersAtlasPropagationCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersAtlasPropagationCollapsibleButton.text = "Atlas Propagation Parameters"
    parametersAtlasPropagationCollapsibleButton.collapsed = True
    self.layout.addWidget(parametersAtlasPropagationCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersAtlasPropagationLayout = qt.QFormLayout(parametersAtlasPropagationCollapsibleButton)

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
    self.setVoxelResolutionLineEditWidget.setText("0.2,0.2,1")
    self.setVoxelResolutionLineEditWidget.setToolTip(
      "Voxel resolution to the image resample function. This is useful to final brain label upsampling which results"
      "in a fine delineation among brain tissues. TIP: Use the typical resolution used in the MRI image acquisition.")
    parametersImageResamplingLayout.addRow("Voxel Resampling Resolution", self.setVoxelResolutionLineEditWidget)

    #
    # Interpolation Functions
    #
    self.setInterpolationFunctionResamplingComboBoxWidget = ctk.ctkComboBox()
    self.setInterpolationFunctionResamplingComboBoxWidget.addItem("bspline")
    self.setInterpolationFunctionResamplingComboBoxWidget.addItem("linear")
    self.setInterpolationFunctionResamplingComboBoxWidget.addItem("nearestNeighbor")
    self.setInterpolationFunctionResamplingComboBoxWidget.addItem("hamming")
    self.setInterpolationFunctionResamplingComboBoxWidget.addItem("cosine")
    self.setInterpolationFunctionResamplingComboBoxWidget.addItem("welch")
    self.setInterpolationFunctionResamplingComboBoxWidget.addItem("lanczos")
    self.setInterpolationFunctionResamplingComboBoxWidget.addItem("blackman")
    self.setInterpolationFunctionResamplingComboBoxWidget.setToolTip(
      "Interpolation functions for resampling step.")
    parametersImageResamplingLayout.addRow("Interpolation ", self.setInterpolationFunctionResamplingComboBoxWidget)

    #
    # Median Filter Parameters Area
    #
    parametersMedianFilterCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersMedianFilterCollapsibleButton.text = "Median Filter Parameters"
    parametersMedianFilterCollapsibleButton.collapsed = True
    self.layout.addWidget(parametersMedianFilterCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersMedianFilterLayout = qt.QFormLayout(parametersMedianFilterCollapsibleButton)

    #
    # Neighborhood Size
    #
    self.setNeighborhoodSizeLineEditWidget = qt.QLineEdit()
    self.setNeighborhoodSizeLineEditWidget.setText("5,5,1")
    self.setNeighborhoodSizeLineEditWidget.setToolTip(
      "Choose the neighborhood applied on the median filter. A large neighborhood will provide a smoother version of the"
      "brain label, however minor details may vanish. TIP: It is commonly used a size of 2x to 5x of the smallest voxel"
      "size, e.g. if the input image has 0.2x0.2x1.0 mm3 (1x1x1 voxels), than the neighborhood size would be [2,2,1].")
    parametersMedianFilterLayout.addRow("Neighborhood Size ", self.setNeighborhoodSizeLineEditWidget)

    #
    # Iterations
    #
    self.setMedianIterationsWidget = qt.QSpinBox()
    self.setMedianIterationsWidget.setMinimum(1)
    self.setMedianIterationsWidget.setMaximum(20)
    self.setMedianIterationsWidget.setValue(2)
    self.setMedianIterationsWidget.setToolTip(
      "Set how many median filtering will be applied. The higher it is, the stronger will be the label smoothing.")
    parametersMedianFilterLayout.addRow("Interations ", self.setMedianIterationsWidget)


    #
    # Apply Brain Volume Refinement
    #
    self.setApplyBrainVolumeRefinementBooleanWidget = ctk.ctkCheckBox()
    self.setApplyBrainVolumeRefinementBooleanWidget.setChecked(True)
    self.setApplyBrainVolumeRefinementBooleanWidget.setToolTip(
      "Check this if you want to refine the brain volume using the atlas brain mask. This step can be helpful if the "
      "brain extraction procedure left some non-brain tissues in the input image.")
    parametersAtlasPropagationLayout.addRow("Apply Brain Volume Refinement", self.setApplyBrainVolumeRefinementBooleanWidget)

    #
    # Brain Atlas
    #
    self.setBrainAtlasComboBoxWidget = ctk.ctkComboBox()
    self.setBrainAtlasComboBoxWidget.addItem("NEO2012") # TODO Ver se usa tambem outro template (2015, http://brain-development.org/brain-atlases/multi-structural-neonatal-brain-atlas/)
    # self.setBrainAtlasComboBoxWidget.addItem("NEO2015") # TODO Novo brain atlas com o mesmo padrao do NEO2012... tem mais detalhes de segmentacao
    self.setBrainAtlasComboBoxWidget.addItem("FET2012") # TODO Preparar cerebellum e brainstem.
    # self.setBrainAtlasComboBoxWidget.addItem("PED2008") # TODO PED2008 precisa separar todas as areas... ver se realmente precisa para agora ou deixa para UPGRADE
    self.setBrainAtlasComboBoxWidget.setToolTip(
      "Choose the most suitable brain atlas for the input image. A list of available atlas are given, however only the "
      "binary labels are considered. These brain atlases will mainly help to segment the cerebellum, brainstem and deep"
      "gray matter. Available atlases: NEO2012 (Neonatal), FET2012 (Fetal) and PED2008 (Pediatric).")
    parametersAtlasPropagationLayout.addRow("Brain Atlas ", self.setBrainAtlasComboBoxWidget)

    #
    # Subject Age
    #
    self.setSubjectAgeIntegerWidget = ctk.ctkSliderWidget()
    self.setSubjectAgeIntegerWidget.maximum = 44
    self.setSubjectAgeIntegerWidget.minimum = 23
    self.setSubjectAgeIntegerWidget.value = 26
    self.setSubjectAgeIntegerWidget.singleStep = 1
    self.setSubjectAgeIntegerWidget.setToolTip("Select the subject's age in weeks. This is only used for neonatal and fetal brain atlases. "
                                               "NOTE: Each atlas has its own age range, with NEO2012=27-44 and FET2012=23-37 weeks, respectivelly."
                                               "If you choose an age below (above), the lower (higher) age will be chosen.")
    parametersAtlasPropagationLayout.addRow("Age ", self.setSubjectAgeIntegerWidget)

    #
    # Registration Algorithm
    #
    self.groupBoxRadioButtons = qt.QGroupBox("Registration Algorithm")
    RadioButtonLayout = qt.QFormLayout()
    self.groupBoxRadioButtons.setLayout(RadioButtonLayout)
    self.setRadioBRAINS = qt.QRadioButton('BRAINSFit')
    self.setRadioBRAINS.setToolTip("Use the Slicer built-in BRAINSFit algorithm (General Registration).")
    self.setRadioANTs = qt.QRadioButton('ANTs')
    self.setRadioANTs.setToolTip("Use the ANTs SyN diffeomorphic algorithm (recommended). If the ANTs tools are not installed in the user's machine, than this option will be not available.")
    if os.environ.get('ANTSPATH'):
      self.setRadioANTs.setChecked(True)
    else:
      self.setRadioBRAINS.setChecked(True)
      self.setRadioANTs.setDisabled(True)
    RadioButtonLayout.addRow(self.setRadioBRAINS)
    RadioButtonLayout.addRow(self.setRadioANTs)

    parametersAtlasPropagationLayout.addRow(self.groupBoxRadioButtons)

    #
    # ANTs Parameters
    #
    self.groupANTsParametersBoxButtons = qt.QGroupBox("ANTs Parameters")
    ANTsParametersLayout = qt.QFormLayout()
    self.groupANTsParametersBoxButtons.setLayout(ANTsParametersLayout)

    #
    # Use Quick Registration
    #
    self.setUseANTSQuickBooleanWidget = ctk.ctkCheckBox()
    self.setUseANTSQuickBooleanWidget.setChecked(False)
    self.setUseANTSQuickBooleanWidget.setToolTip(
      "Check this if you want to use the antsRegistrationSyNQuick.sh script. This will considerably reduce the "
      "total time required in the registration step.")

    #
    # Number of Cores
    #
    self.setNumberOfCoresWidget = ctk.ctkSliderWidget()
    self.setNumberOfCoresWidget.singleStep = 1
    self.setNumberOfCoresWidget.minimum = 1
    self.setNumberOfCoresWidget.maximum = multiprocessing.cpu_count()
    self.setNumberOfCoresWidget.value = self.setNumberOfCoresWidget.maximum - 1
    self.setNumberOfCoresWidget.setToolTip(
      "Set the number of CPU's used in the registration process. In order to prevent the SO crash, it is advisable to use N - 1 (N = Total number of cores available).")

    #
    # Radius for correlation calculation
    #
    self.setCorrelationRadiusWidget = qt.QSpinBox()
    self.setCorrelationRadiusWidget.setMinimum(1)
    self.setCorrelationRadiusWidget.setMaximum(30)
    self.setCorrelationRadiusWidget.setValue(4)
    self.setCorrelationRadiusWidget.setToolTip(
      "Set the radius for cross correlation metric used in the SyN registration. Units given in number of voxels.")

    ANTsParametersLayout.addRow("Use Quick Registration ", self.setUseANTSQuickBooleanWidget)
    ANTsParametersLayout.addRow("Number Of Cores ", self.setNumberOfCoresWidget)
    ANTsParametersLayout.addRow("Radius ", self.setCorrelationRadiusWidget)

    parametersAtlasPropagationLayout.addRow(self.groupANTsParametersBoxButtons)

    #
    # BRAINSFit Parameters
    #
    self.groupBRAINFitParametersBoxButtons = qt.QGroupBox("BRAINSFit Parameters")
    BRAINSFitParametersLayout = qt.QFormLayout()
    self.groupBRAINFitParametersBoxButtons.setLayout(BRAINSFitParametersLayout)

    #
    # Percentage Sampling Area
    #
    self.setPercSamplingQWidget = qt.QDoubleSpinBox()
    self.setPercSamplingQWidget.setDecimals(4)
    self.setPercSamplingQWidget.setMaximum(0.1)
    self.setPercSamplingQWidget.setMinimum(0.0001)
    self.setPercSamplingQWidget.setSingleStep(0.0001)
    self.setPercSamplingQWidget.setValue(0.05)
    self.setPercSamplingQWidget.setToolTip("Percentage of voxel used in registration.")

    #
    # BSpline Grid
    #
    self.setBSplineGridWidget = qt.QLineEdit()
    self.setBSplineGridWidget.setText('14,14,10')
    self.setBSplineGridWidget.setToolTip("Set the BSpline grid for non linear structural adjustments.")

    #
    # Initiation Method Area
    #
    self.setInitiationRegistrationBooleanWidget = ctk.ctkComboBox()
    self.setInitiationRegistrationBooleanWidget.addItem("Off")
    self.setInitiationRegistrationBooleanWidget.addItem("useMomentsAlign")
    self.setInitiationRegistrationBooleanWidget.addItem("useCenterOfHeadAlign")
    self.setInitiationRegistrationBooleanWidget.addItem("useGeometryAlign")
    self.setInitiationRegistrationBooleanWidget.setToolTip(
      "Initialization method used for the MNI152 registration.")

    #
    # Cost Metric
    #
    self.setRegistrationCostMetricWidget = ctk.ctkComboBox()
    self.setRegistrationCostMetricWidget.addItem("NC")
    self.setRegistrationCostMetricWidget.addItem("MMI")
    self.setRegistrationCostMetricWidget.addItem("MSE")
    self.setRegistrationCostMetricWidget.setToolTip(
      "The cost metric to be used during fitting. Defaults to NC. Options are MMI (Mattes Mutual Information), MSE (Mean Square Error), NC"
      " (Normalized Correlation).")

    #
    # Interpolation
    #
    self.setInterpolationFunctionRegistrationComboBoxWidget = ctk.ctkComboBox()
    self.setInterpolationFunctionRegistrationComboBoxWidget.addItem("BSpline")
    self.setInterpolationFunctionRegistrationComboBoxWidget.addItem("Linear")
    self.setInterpolationFunctionRegistrationComboBoxWidget.addItem("NearestNeighbor")
    self.setInterpolationFunctionRegistrationComboBoxWidget.addItem("WindowedSinc")
    self.setInterpolationFunctionRegistrationComboBoxWidget.addItem("Cosine")
    self.setInterpolationFunctionRegistrationComboBoxWidget.addItem("Welch")
    self.setInterpolationFunctionRegistrationComboBoxWidget.addItem("Lanczos")
    self.setInterpolationFunctionRegistrationComboBoxWidget.addItem("Blackman")
    self.setInterpolationFunctionRegistrationComboBoxWidget.addItem("Hamming")
    self.setInterpolationFunctionRegistrationComboBoxWidget.setToolTip(
      "Interpolation functions for registration step.")

    BRAINSFitParametersLayout.addRow("Percentage Of Samples ", self.setPercSamplingQWidget)
    BRAINSFitParametersLayout.addRow("Spline Grid ", self.setBSplineGridWidget)
    BRAINSFitParametersLayout.addRow("Initiation Method ", self.setInitiationRegistrationBooleanWidget)
    BRAINSFitParametersLayout.addRow("Cost Metric ", self.setRegistrationCostMetricWidget)
    BRAINSFitParametersLayout.addRow("Interpolation ", self.setInterpolationFunctionRegistrationComboBoxWidget)

    parametersAtlasPropagationLayout.addRow(self.groupBRAINFitParametersBoxButtons)


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
    self.brainTissuesSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSelect()

  def cleanup(self):
    pass

  def onSelect(self):
    self.applyButton.enabled = self.inputSelector.currentNode() and self.brainTissuesSelector.currentNode()

  def onApplyButton(self):
    logic = BabyBrainSegmentationLogic()
    modality = self.setImageModalityBooleanWidget.currentText
    estimateBasalGanglia = self.setUseBasalGangliaEstimatorBooleanWidget.isChecked()
    # debugMode = self.setUseDebugModeBooleanWidget.isChecked()
    applyBrainVolumeRefinements = self.setApplyBrainVolumeRefinementBooleanWidget.isChecked()
    brainAtlas = self.setBrainAtlasComboBoxWidget.currentText
    age = self.setSubjectAgeIntegerWidget.value
    if self.setRadioBRAINS.isChecked():
      registrationAlgorithm = self.setRadioBRAINS.text
    else:
      registrationAlgorithm = self.setRadioANTs.text
    useQuickRegistration = self.setUseANTSQuickBooleanWidget.isChecked()
    numOfCores = self.setNumberOfCoresWidget.value
    correlationRadius = self.setCorrelationRadiusWidget.value
    sampling = self.setPercSamplingQWidget.value
    splineGrid = self.setBSplineGridWidget.text
    initMethod = self.setInitiationRegistrationBooleanWidget.currentText
    costMetric = self.setRegistrationCostMetricWidget.currentText
    interpolationRegistration = self.setInterpolationFunctionRegistrationComboBoxWidget.currentText
    voxelResampling = self.setVoxelResolutionLineEditWidget.text
    interpolationResampling = self.setInterpolationFunctionResamplingComboBoxWidget.currentText
    neighborSize = self.setNeighborhoodSizeLineEditWidget.text
    interations = self.setMedianIterationsWidget.value
    logic.run(self.inputSelector.currentNode()
              , self.brainTissuesSelector.currentNode()
              , modality
              , estimateBasalGanglia
              # , debugMode
              , applyBrainVolumeRefinements
              , brainAtlas
              , age
              , registrationAlgorithm
              , useQuickRegistration
              , numOfCores
              , correlationRadius
              , sampling
              , splineGrid
              , initMethod
              , costMetric
              , interpolationRegistration
              , voxelResampling
              , interpolationResampling
              , neighborSize
              , interations)

#
# BabyBrainSegmentationLogic
#

class BabyBrainSegmentationLogic(ScriptedLoadableModuleLogic):
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
          , modality
          , estimateBasalGanglia
          # , debugMode
          , applyBrainVolumeRefinements
          , brainAtlas
          , age
          , registrationAlgorithm
          , useQuickRegistration
          , numOfCores
          , correlationRadius
          , sampling
          , splineGrid
          , initMethod
          , costMetric
          , interpolationRegistration
          , voxelResampling
          , interpolationResampling
          , neighborSize
          , interations):
    """
    Run the actual algorithm
    """

    if not self.isValidInputOutputData(inputVolume, outputVolume):
      slicer.util.errorDisplay('Input volume is the same as output volume. Choose a different output volume.')
      return False

    if brainAtlas == "FET2012" and modality == "T1":
        slicer.util.errorDisplay('The FET2012 atlas does not have a T1 brain template. Choose a different image modality or a different brain template.')
        return False

    if brainAtlas == "FET2012" and estimateBasalGanglia:
        slicer.util.errorDisplay('The FET2012 atlas does not have an option to estimate deep gray matter structures. Choose a different brain template or leave unchecked the Estimate Basal Ganglia Structures option .')
        return False

    modulePath = os.path.dirname(slicer.modules.babybrainsegmentation.path)
    if platform.system() is "Windows":
      databasePath = modulePath + "\\Resources\\atlases"
    else:
      databasePath = modulePath + "/Resources/atlases"

    # TODO Testes que ja foram feitos:
    # 1) Usar transformadas faz com que o atlas propagation pule corretamente
    # 2)

    #  TODO Fazer um checkup pra ver se este color table ja foi carregado no slicer...assim nao carrega duas vezes
    # Loading BabyBrain color table
    if platform.system() is "Windows":
      slicer.util.loadColorTable(modulePath + "\\Resources\\2018_USP_BabyBrain.ctbl")
    else:
      slicer.util.loadColorTable(modulePath + "/Resources/2018_USP_BabyBrain.ctbl")

    logging.info('Processing started')

    # Checking if the age is found in the chosen brain template
    setAge = age;
    if brainAtlas == "NEO2012":
        if age < 27:
            setAge = 27
        elif age > 43:
            setAge = 43
    elif brainAtlas == "FET2012":
        if age < 23:
            setAge = 23
        elif age > 37:
            setAge = 37

    # This will check is already exist the registration transformation in the Slicer scene, using the generic names of it.
    # If those exist, the most recent one is adopted as the correct template transformations to the input data.
    # This strategy is useful if the user used the BabyBrainPrepation module, because the registration transformation
    # generated there are recycled here.
    jumpedAtlasPropagation = False
    regAffine=slicer.util.getNodes('BabyBrain_regMNI2Native_0GenericAffine')
    regWarp=slicer.util.getNodes('BabyBrain_regMNI2Native_1Warp')
    if regAffine.get('BabyBrain_regMNI2Native_0GenericAffine') is None and regWarp.get('BabyBrain_regMNI2Native_1Warp') is None:
        jumpedAtlasPropagation = True
        ######################################################################################
        # Step  - Label propagation using brain atlas.
        ######################################################################################
        # Reading brain template
        if platform.system() is "Windows":
          readingParameters = {}
          readingParameters['name'] = "brain_template"
          readingParameters['center'] = True
          readingParameters['show'] = False
          (readSuccess, brainAtlasNode) = slicer.util.loadVolume(databasePath +
                                                                 "\\" + brainAtlas +
                                                                 "\\templates\\template_" + modality + "_" + str(int(setAge)) +".nii.gz", readingParameters, True)
        else:
          readingParameters = {}
          readingParameters['name'] = "brain_template"
          readingParameters['center'] = True
          readingParameters['show'] = False
          (readSuccess, brainAtlasNode) = slicer.util.loadVolume(databasePath +
                                                                 "/" + brainAtlas +
                                                                 "/templates/template_" + modality + "_" + str(int(setAge)) + ".nii.gz", readingParameters, True)

        ######################################################################################
        # Step  - Atlas propagation - linear and elastic transformations
        ######################################################################################
        # Image registration with atlas - ANTs or BRAINSfit
        # Creating linear transform node
        regMNI2NativeLinearTransform = slicer.vtkMRMLLinearTransformNode()
        regMNI2NativeLinearTransform.SetName("BabyBrain_regMNI2Native_0GenericAffine")
        slicer.mrmlScene.AddNode(regMNI2NativeLinearTransform)

        regMNI2NativeBSplineTransform = slicer.vtkMRMLBSplineTransformNode()
        regMNI2NativeBSplineTransform.SetName("BabyBrain_regMNI2Native_1Warp")
        slicer.mrmlScene.AddNode(regMNI2NativeBSplineTransform)


        self.atlasPropagation(registrationAlgorithm
                              , inputVolume
                              , brainAtlasNode
                              , regMNI2NativeLinearTransform
                              , regMNI2NativeBSplineTransform
                              , interpolationRegistration
                              , sampling
                              , splineGrid
                              , initMethod
                              , numOfCores
                              , useQuickRegistration)

        if registrationAlgorithm == "ANTs":
          slicer.mrmlScene.RemoveNode(regMNI2NativeLinearTransform)
          slicer.mrmlScene.RemoveNode(regMNI2NativeBSplineTransform)

    # Allocating the default path in order to save temporary files into the hard drive
    home = expanduser("~") # TODO colocar home para windows tambem ... expanduser(%USERPROFILE%) algo assim...
    tmpFolder=home
    # Creating temporary folder in home directory
    if registrationAlgorithm == "ANTs":
        tmpFolder = home + "/tmpANTsBabyBrainSegmentation"
    else:
        tmpFolder = home + "/tmpBRAINSFitBabyBrainSegmentation" # TODO Colocar uma pasta temporaria para o BRAINSFit tambem!! fazer com isso seja o output geral da extensao, que o usuario pode querer salvar se quiser!!! (Debug mode foi parecido com isso!!)




    slicer.util.showStatusMessage("Atlas propagation is finished...")
    ######################################################################################
    # Step  - Resampling the input volume
    ######################################################################################
    tmpResampledInputNode = slicer.vtkMRMLScalarVolumeNode()
    tmpResampledInputNode.SetName("resampled_input")
    slicer.mrmlScene.AddNode(tmpResampledInputNode)
    self.imageResamplingResolution(inputVolume
                                   , tmpResampledInputNode
                                   , voxelResampling
                                   , interpolationResampling)
    slicer.util.showStatusMessage("Image resampling to voxel resolution of [" + str(voxelResampling) + "] is finished...")

    # Applying a brain volume refinement using atlas brain mask (if required) # TODO Melhorar este brain volume correction... pensar em como reduzir os erros do BET aqui...criar um cli especfico para isto
    if applyBrainVolumeRefinements:
      # Reading brain volume mask from atlas
      if platform.system() is "Windows":
        readingParameters = {}
        readingParameters['name'] = "brain_template_mask"
        readingParameters['center'] = True
        readingParameters['show'] = False
        readingParameters['labelmap'] = True
        (readSuccess, brainMaskNode) = slicer.util.loadVolume(databasePath +
                                                               "\\" + brainAtlas +
                                                               "\\brainmask\\brainmask_" + str(int(setAge)) + ".nii.gz", readingParameters,
                                                               True)
      else:
        readingParameters = {}
        readingParameters['name'] = "brain_template_mask"
        readingParameters['center'] = True
        readingParameters['show'] = False
        readingParameters['labelmap'] = True
        (readSuccess, brainMaskNode) = slicer.util.loadVolume(databasePath +
                                                               "/" + brainAtlas +
                                                               "/brainmask/brainmask_" + str(int(setAge)) + ".nii.gz", readingParameters,
                                                               True)

      # Transforming brain mask to native space
      tmpBrainMask = slicer.vtkMRMLLabelMapVolumeNode()
      tmpBrainMask.SetName("brainmask")
      slicer.mrmlScene.AddNode(tmpBrainMask)
      self.applyRegistrationTransforms(registrationAlgorithm
                                       , brainMaskNode
                                       , tmpResampledInputNode
                                       , tmpBrainMask
                                       , slicer.util.getNode("BabyBrain_regMNI2Native_0GenericAffine")
                                       , slicer.util.getNode("BabyBrain_regMNI2Native_1Warp")
                                       , True)

      # Applying brain volume correction
      params = {}
      params['InputVolume'] = tmpResampledInputNode.GetID()
      params['MaskVolume'] = tmpBrainMask.GetID()
      params['OutputVolume'] = tmpResampledInputNode.GetID()
      params['Label'] = 1

      slicer.cli.run(slicer.modules.maskscalarvolume, None, params, wait_for_completion=True)
      slicer.util.showStatusMessage("Brain volume refinements is finished...")

      # Removing tmp files
      # slicer.mrmlScene.RemoveNode(tmpBrainMask)
      # slicer.mrmlScene.RemoveNode(brainMaskNode)

    ######################################################################################
    # Step  - Remove Cerebellum and brainstem from the input volume
    ######################################################################################
    # Reading the cerebellum region from atlas and adjusting to native space
    tmpCerebellumMask = slicer.vtkMRMLLabelMapVolumeNode()
    tmpCerebellumMask.SetName("cerebellum_mask")
    slicer.mrmlScene.AddNode(tmpCerebellumMask)
    # Reading cerebellum volume mask from atlas
    if platform.system() is "Windows":
      readingParameters = {}
      readingParameters['name'] = "cerebellum_template_mask"
      readingParameters['center'] = True
      readingParameters['show'] = False
      readingParameters['labelmap'] = True
      (readSuccess, cerebellumMaskNode) = slicer.util.loadVolume(databasePath +
                                                            "\\" + brainAtlas +
                                                            "\\cerebellum\\cerebellum_" + str(int(setAge)) + ".nii.gz", readingParameters,
                                                            True)
    else:
      readingParameters = {}
      readingParameters['name'] = "cerebellum_template_mask"
      readingParameters['center'] = True
      readingParameters['show'] = False
      readingParameters['labelmap'] = True
      (readSuccess, cerebellumMaskNode) = slicer.util.loadVolume(databasePath +
                                                            "/" + brainAtlas +
                                                            "/cerebellum/cerebellum_" + str(int(setAge)) + ".nii.gz", readingParameters,
                                                            True)

    self.applyRegistrationTransforms(registrationAlgorithm
                                     , cerebellumMaskNode
                                     , tmpResampledInputNode
                                     , tmpCerebellumMask
                                     , slicer.util.getNode("BabyBrain_regMNI2Native_0GenericAffine")
                                     , slicer.util.getNode("BabyBrain_regMNI2Native_1Warp")
                                     , True)

    # Reading the cerebellum prior probability and adjusting to native space
    tmpCerebellumPriors = slicer.vtkMRMLScalarVolumeNode()
    tmpCerebellumPriors.SetName("cerebellum_priors")
    slicer.mrmlScene.AddNode(tmpCerebellumPriors)
    # Reading cerebellum volume mask from atlas
    if platform.system() is "Windows":
        readingParameters = {}
        readingParameters['name'] = "cerebellum_template_priors"
        readingParameters['center'] = True
        readingParameters['show'] = False
        readingParameters['labelmap'] = False
        (readSuccess, cerebellumPriorsNode) = slicer.util.loadVolume(databasePath +
                                                                   "\\" + brainAtlas +
                                                                   "\\cerebellum\\cerebellum_prob_" + str(int(setAge)) + ".nii.gz", readingParameters,
                                                                   True)
    else:
        readingParameters = {}
        readingParameters['name'] = "cerebellum_template_priors"
        readingParameters['center'] = True
        readingParameters['show'] = False
        readingParameters['labelmap'] = False
        (readSuccess, cerebellumPriorsNode) = slicer.util.loadVolume(databasePath +
                                                                   "/" + brainAtlas +
                                                                   "/cerebellum/cerebellum_prob_" + str(int(setAge)) + ".nii.gz", readingParameters,
                                                                   True)

    self.applyRegistrationTransforms(registrationAlgorithm
                                     , cerebellumPriorsNode
                                     , tmpResampledInputNode
                                     , tmpCerebellumPriors
                                     , slicer.util.getNode("BabyBrain_regMNI2Native_0GenericAffine")
                                     , slicer.util.getNode("BabyBrain_regMNI2Native_1Warp")
                                     , False)

    # Saving the cerebellum priors image
    slicer.util.saveNode(tmpCerebellumPriors, tmpFolder + '/cerebellum_priors.nii.gz')


    # Reading the brainstem region from atlas and adjusting to native space
    tmpBrainstemMask = slicer.vtkMRMLLabelMapVolumeNode()
    tmpBrainstemMask.SetName("brainstem_mask")
    slicer.mrmlScene.AddNode(tmpBrainstemMask)
    # Reading brainstem volume mask from atlas
    if platform.system() is "Windows":
      readingParameters = {}
      readingParameters['name'] = "brainstem_template_mask"
      readingParameters['center'] = True
      readingParameters['show'] = False
      readingParameters['labelmap'] = True
      (readSuccess, brainstemMaskNode) = slicer.util.loadVolume(databasePath +
                                                                 "\\" + brainAtlas +
                                                                 "\\brainstem\\brainstem_" + str(int(setAge)) + ".nii.gz", readingParameters,
                                                                 True)
    else:
      readingParameters = {}
      readingParameters['name'] = "brainstem_template_mask"
      readingParameters['center'] = True
      readingParameters['show'] = False
      readingParameters['labelmap'] = True
      (readSuccess, brainstemMaskNode) = slicer.util.loadVolume(databasePath +
                                                                 "/" + brainAtlas +
                                                                 "/brainstem/brainstem_" + str(int(setAge)) + ".nii.gz", readingParameters,
                                                                 True)

    self.applyRegistrationTransforms(registrationAlgorithm
                                     , brainstemMaskNode
                                     , tmpResampledInputNode
                                     , tmpBrainstemMask
                                     , slicer.util.getNode("BabyBrain_regMNI2Native_0GenericAffine")
                                     , slicer.util.getNode("BabyBrain_regMNI2Native_1Warp")
                                     , True)

    # Reading the brainstem priors and adjusting to native space
    tmpBrainstemPriors = slicer.vtkMRMLScalarVolumeNode()
    tmpBrainstemPriors.SetName("brainstem_priors")
    slicer.mrmlScene.AddNode(tmpBrainstemPriors)
    # Reading brainstem volume mask from atlas
    if platform.system() is "Windows":
        readingParameters = {}
        readingParameters['name'] = "brainstem_template_priors"
        readingParameters['center'] = True
        readingParameters['show'] = False
        readingParameters['labelmap'] = False
        (readSuccess, brainstemPriorsNode) = slicer.util.loadVolume(databasePath +
                                                                  "\\" + brainAtlas +
                                                                  "\\brainstem\\brainstem_prob_" + str(int(setAge)) + ".nii.gz", readingParameters,
                                                                  True)
    else:
        readingParameters = {}
        readingParameters['name'] = "brainstem_template_priors"
        readingParameters['center'] = True
        readingParameters['show'] = False
        readingParameters['labelmap'] = False
        (readSuccess, brainstemPriorsNode) = slicer.util.loadVolume(databasePath +
                                                                  "/" + brainAtlas +
                                                                  "/brainstem/brainstem_prob_" + str(int(setAge)) + ".nii.gz", readingParameters,
                                                                  True)

    self.applyRegistrationTransforms(registrationAlgorithm
                                     , brainstemPriorsNode
                                     , tmpResampledInputNode
                                     , tmpBrainstemPriors
                                     , slicer.util.getNode("BabyBrain_regMNI2Native_0GenericAffine")
                                     , slicer.util.getNode("BabyBrain_regMNI2Native_1Warp")
                                     , False)

    # Saving the brainstem priors image
    slicer.util.saveNode(tmpBrainstemPriors, tmpFolder + '/brainstem_priors.nii.gz')


    # Removing brainstem and cerebellum from the input image
    inputImage = sitkUtils.PullVolumeFromSlicer(tmpResampledInputNode)
    brainstemMaskImage = sitkUtils.PullVolumeFromSlicer(tmpBrainstemMask)
    brainstemMaskImage = sitk.Cast(brainstemMaskImage, inputImage.GetPixelID()) # This is necessary since the MaskNegated filter requires the same pixel type in both images.
    cerebellumMaskImage = sitkUtils.PullVolumeFromSlicer(tmpCerebellumMask)
    cerebellumMaskImage = sitk.Cast(cerebellumMaskImage, inputImage.GetPixelID())  # This is necessary since the MaskNegated filter requires the same pixel type in both images.
    filter = sitk.MaskNegatedImageFilter()
    output_brainOnly_Image = filter.Execute(inputImage,brainstemMaskImage)
    output_brainOnly_Image = filter.Execute(output_brainOnly_Image, cerebellumMaskImage)

    tmpBrainOnlyNode = slicer.vtkMRMLScalarVolumeNode()
    tmpBrainOnlyNode.SetName("resampled_input_brainOnly")
    slicer.mrmlScene.AddNode(tmpBrainOnlyNode)
    sitkUtils.PushVolumeToSlicer(output_brainOnly_Image, tmpBrainOnlyNode) # TODO Ver como carregar estas labels sem atualizar a scene... nao fica legal ver essas labels intermediarias durante o processamento

    slicer.util.showStatusMessage("Brainstem and cerebellum removal is finished...")

    ######################################################################################
    # Step  - Segmenting the cerebellum and brainstem # TODO Pode colocar as priors do atlas aqui...so fazer o resampling com interpolacao linear...atualizar o Bayesian tissue classifier module para poder entrar com o prior
    ######################################################################################
    # Cerebellum
    tmpCerebellumOnlyVolumeNode = slicer.vtkMRMLScalarVolumeNode()
    tmpCerebellumOnlyVolumeNode.SetName("cerebellum_volume")
    slicer.mrmlScene.AddNode(tmpCerebellumOnlyVolumeNode)
    # Extracting the cerebellum from the input image
    params = {}
    params['InputVolume'] = tmpResampledInputNode.GetID()
    params['MaskVolume'] = tmpCerebellumMask.GetID()
    params['OutputVolume'] = tmpCerebellumOnlyVolumeNode.GetID()
    params['Label'] = 1

    slicer.cli.run(slicer.modules.maskscalarvolume, None, params, wait_for_completion=True)

    cerebellumOnlyLabelMask = slicer.vtkMRMLLabelMapVolumeNode()
    cerebellumOnlyLabelMask.SetName("cerebellum_tissue_mask")
    slicer.mrmlScene.AddNode(cerebellumOnlyLabelMask)
    definitions = [5,3]
    if modality == 'T1':
      definitions = [3,5]

    path2CerebellumPriors=tmpFolder + '/cerebellum_priors.nii.gz' # TODO ERRO DE SEGMENTACAO AQUI PARA FRENTE...NAO SAI COM PRIORS CORRETAMENTE!!
    self.segmentingTissues(tmpCerebellumOnlyVolumeNode
                           , modality
                           , 3
                           , cerebellumOnlyLabelMask
                           , inputPriorsFile=path2CerebellumPriors
                           , definitions=definitions ) # Cerebellum and CSF

    # Brainstem
    tmpBrainstemOnlyVolumeNode = slicer.vtkMRMLScalarVolumeNode()
    tmpBrainstemOnlyVolumeNode.SetName("brainstem_volume")
    slicer.mrmlScene.AddNode(tmpBrainstemOnlyVolumeNode)
    # Extracting the brainstem from the input image
    params = {}
    params['InputVolume'] = tmpResampledInputNode.GetID()
    params['MaskVolume'] = tmpBrainstemMask.GetID()
    params['OutputVolume'] = tmpBrainstemOnlyVolumeNode.GetID()
    params['Label'] = 1

    slicer.cli.run(slicer.modules.maskscalarvolume, None, params, wait_for_completion=True)

    brainstemOnlyLabelMask = slicer.vtkMRMLLabelMapVolumeNode()
    brainstemOnlyLabelMask.SetName("brainstem_tissue_mask")
    slicer.mrmlScene.AddNode(brainstemOnlyLabelMask)
    definitions = [4,3]
    if modality == 'T1':
      definitions = [3,4]

    path2BrainstemPriors = tmpFolder + '/brainstem_priors.nii.gz'
    self.segmentingTissues(tmpBrainstemOnlyVolumeNode
                           , modality
                           , 3
                           , brainstemOnlyLabelMask
                           , inputPriorsFile=path2BrainstemPriors
                           , definitions=definitions) # Brainstem and CSF

    # Merging both labels into one
    brainstemPlusCerebellumLabelMask = slicer.vtkMRMLLabelMapVolumeNode()
    brainstemPlusCerebellumLabelMask.SetName("brainstem+cerebellum_tissue_mask")
    slicer.mrmlScene.AddNode(brainstemPlusCerebellumLabelMask)
    self.combineLabels(brainstemOnlyLabelMask, cerebellumOnlyLabelMask, brainstemPlusCerebellumLabelMask)

    slicer.util.showStatusMessage("Brainstem and cerebellum segmentation is finished...")

    # ######################################################################################
    # # Step  - Brain tissue segmentation (only brain)
    # ######################################################################################
    # brainOnlyLabelMask = slicer.vtkMRMLLabelMapVolumeNode()
    # brainOnlyLabelMask.SetName("brainOnly_tissues_mask")
    # slicer.mrmlScene.AddNode(brainOnlyLabelMask)
    #
    # definitions = [1,2,3]
    # if modality == 'T1':
    #   definitions = [3,2,1]
    #
    # self.segmentingTissues(tmpBrainOnlyNode
    #                        , modality
    #                        , 4
    #                        , brainOnlyLabelMask
    #                        , definitions=definitions)
    # slicer.util.showStatusMessage("Whole brain (first) segmentation is finished...")
    #
    # ######################################################################################
    # # Step  - Smoothing brain segmentation by median filtering
    # ######################################################################################
    # for i in range(0, interations):
    #   params = {}
    #   params['neighborhood'] = neighborSize
    #   params['inputVolume'] = brainOnlyLabelMask
    #   params['outputVolume'] = brainOnlyLabelMask
    #
    #   slicer.cli.run(slicer.modules.medianimagefilter, None, params, wait_for_completion=True)
    #   slicer.util.showStatusMessage("Brain segmentation smoothing is finished...")
    #

    #
    # ######################################################################################
    # # Step  - Split brain hemispheres (only adjust the mask on the native space)
    # ######################################################################################
    # # Reading brain hemispheres mask from atlas
    # if platform.system() is "Windows":
    #   readingParameters = {}
    #   readingParameters['name'] = "hemisphere_template_mask"
    #   readingParameters['center'] = True
    #   readingParameters['show'] = False
    #   readingParameters['labelmap'] = True
    #   (readSuccess, hemispheresMaskNode) = slicer.util.loadVolume(databasePath +
    #                                                              "\\" + brainAtlas +
    #                                                              "\\hemispheres.nii.gz", readingParameters,
    #                                                              True)
    # else:
    #   readingParameters = {}
    #   readingParameters['name'] = "hemisphere_template_mask"
    #   readingParameters['center'] = True
    #   readingParameters['show'] = False
    #   readingParameters['labelmap'] = True
    #   (readSuccess, hemispheresMaskNode) = slicer.util.loadVolume(databasePath +
    #                                                              "/" + brainAtlas +
    #                                                              "/hemispheres.nii.gz", readingParameters,
    #                                                              True)
    #
    # # Transforming the hemispheres mask to native space
    # tmpHemispheresLabelMask = slicer.vtkMRMLLabelMapVolumeNode()
    # tmpHemispheresLabelMask.SetName("tmpHemispheres_mask")
    # slicer.mrmlScene.AddNode(tmpHemispheresLabelMask)
    # self.applyRegistrationTransforms(registrationAlgorithm
    #                                  , hemispheresMaskNode
    #                                  , tmpResampledInputNode
    #                                  , tmpHemispheresLabelMask
    #                                  , slicer.util.getNode("BabyBrain_regMNI2Native_0GenericAffine")
    #                                  , slicer.util.getNode("BabyBrain_regMNI2Native_1Warp")
    #                                  , True)
    # slicer.util.showStatusMessage("Brain hemispheres space definition is finished...")
    #
    # ######################################################################################
    # # Step  - Setting up the brain hemispheres value (white matter, gray matter and ventricules)
    # ######################################################################################
    # # Reading brain hemispheres mask from atlas
    # if platform.system() is "Windows":
    #     readingParameters = {}
    #     readingParameters['name'] = "ventricules_template_mask"
    #     readingParameters['center'] = True
    #     readingParameters['show'] = False
    #     readingParameters['labelmap'] = True
    #     (readSuccess, venctriculesMaskNode) = slicer.util.loadVolume(databasePath +
    #                                                                 "\\" + brainAtlas +
    #                                                                 "\\ventricules\\ventricules_" + str(int(setAge)) + ".nii.gz", readingParameters,
    #                                                                 True)
    # else:
    #     readingParameters = {}
    #     readingParameters['name'] = "ventricules_template_mask"
    #     readingParameters['center'] = True
    #     readingParameters['show'] = False
    #     readingParameters['labelmap'] = True
    #     (readSuccess, venctriculesMaskNode) = slicer.util.loadVolume(databasePath +
    #                                                                 "/" + brainAtlas +
    #                                                                 "/ventricules/ventricules_" + str(int(setAge)) + ".nii.gz", readingParameters,
    #                                                                 True)
    #
    # # Transforming the hemispheres mask to native space
    # tmpVentriculesLabelMask = slicer.vtkMRMLLabelMapVolumeNode()
    # tmpVentriculesLabelMask.SetName("tmpVentricules_mask")
    # slicer.mrmlScene.AddNode(tmpVentriculesLabelMask)
    # self.applyRegistrationTransforms(registrationAlgorithm
    #                                  , venctriculesMaskNode
    #                                  , tmpResampledInputNode
    #                                  , tmpVentriculesLabelMask
    #                                  , slicer.util.getNode("BabyBrain_regMNI2Native_0GenericAffine")
    #                                  , slicer.util.getNode("BabyBrain_regMNI2Native_1Warp")
    #                                  , True)
    #
    # # Adjusting the brain hemispheres labels, taking care to ventricules being set to value 10
    # brainOnlyHemispheresLabelMask = slicer.vtkMRMLLabelMapVolumeNode()
    # brainOnlyHemispheresLabelMask.SetName("brain_hemispheres_tissue_mask")
    # slicer.mrmlScene.AddNode(brainOnlyHemispheresLabelMask)
    # params = {}
    # params['inputLabel'] = brainOnlyLabelMask.GetID()
    # params['splitLabel'] = tmpVentriculesLabelMask.GetID()
    # params['outputLabel'] = brainOnlyHemispheresLabelMask.GetID() # TODO Verificar se ao inves de marcar os ventriculos apenas pelos labels, fazer a segmentacao novamente pode ser melhor. O cerebro pode ja vir filtrado com median filter e os ventriculos sao resegmentados...
    # if modality == "T1":
    #     params['labelSideA'] = 9
    # else:
    #     params['labelSideA'] = 7
    # params['labelSideB'] = 0
    # params['doKeepSomeValues'] = True
    # if modality == 'T1':
    #     params['keepSideA'] = "2,3"
    #     params['keepSideB'] = "1,2,3"
    # else:
    #     params['keepSideA'] = "1,2"
    #     params['keepSideB'] = "1,2,3"
    #
    # slicer.cli.run(slicer.modules.splitlabelvalues, None, params, wait_for_completion=True)
    # slicer.util.showStatusMessage("Brain hemispheres detection and separation is finished...")
    #
    # # Adjusting the brain hemispheres globally, i.e. right-left white, gray and ventricules. CSF is constant regardless laterality
    # params = {}
    # params['inputLabel'] = brainOnlyHemispheresLabelMask.GetID()
    # params['splitLabel'] = tmpHemispheresLabelMask.GetID()
    # params['outputLabel'] = brainOnlyHemispheresLabelMask.GetID()
    # params['labelSideA'] = 20
    # params['labelSideB'] = 40
    # params['doKeepSomeValues'] = True  # This will maintain the CSF intact
    # csf = 3
    # if modality == 'T1':
    #   csf = 1
    # params['keepSideA'] = csf
    # params['keepSideB'] = csf
    #
    # slicer.cli.run(slicer.modules.splitlabelvalues, None, params, wait_for_completion=True)
    #
    # slicer.util.showStatusMessage("Ventricules detection and separation is finished...")
    #
    # ######################################################################################
    # # Step  - Merging brainstem, cerebellum and brain hemispheres all the tissues together
    # ######################################################################################
    # fullBrainSegmentationLabelMask = slicer.vtkMRMLLabelMapVolumeNode()
    # fullBrainSegmentationLabelMask.SetName("full_brain_tissues_mask")
    # slicer.mrmlScene.AddNode(fullBrainSegmentationLabelMask)
    # self.combineLabels(brainOnlyHemispheresLabelMask, brainstemPlusCerebellumLabelMask, outputVolume)
    #
    # slicer.util.showStatusMessage("Brain parcellation is finished...")
    #
    # if estimateBasalGanglia:
    #   ######################################################################################
    #   # Step  - Merging basal ganglia to other labels all the tissues together
    #   ######################################################################################
    #   # Reading brain hemispheres mask from atlas
    #   if platform.system() is "Windows":
    #     readingParameters = {}
    #     readingParameters['name'] = "deepGM_template_mask"
    #     readingParameters['center'] = True
    #     readingParameters['show'] = False
    #     readingParameters['labelmap'] = True
    #     (readSuccess, dgmMaskNode) = slicer.util.loadVolume(databasePath +
    #                                                                 "\\" + brainAtlas +
    #                                                                 "\\dgm\\dgm_" + str(int(setAge)) + ".nii.gz", readingParameters,
    #                                                                 True)
    #   else:
    #     readingParameters = {}
    #     readingParameters['name'] = "deepGM_template_mask"
    #     readingParameters['center'] = True
    #     readingParameters['show'] = False
    #     readingParameters['labelmap'] = True
    #     (readSuccess, dgmMaskNode) = slicer.util.loadVolume(databasePath +
    #                                                                 "/" + brainAtlas +
    #                                                                 "/dgm/dgm_" + str(int(setAge)) + ".nii.gz", readingParameters,
    #                                                                 True)
    #
    #   # Transforming the basal ganglia mask to native space
    #   tmpDeepGMLabelMask = slicer.vtkMRMLLabelMapVolumeNode()
    #   tmpDeepGMLabelMask.SetName("tmpdeepGM_mask")
    #   slicer.mrmlScene.AddNode(tmpDeepGMLabelMask)
    #   self.applyRegistrationTransforms(registrationAlgorithm
    #                                    , dgmMaskNode
    #                                    , tmpResampledInputNode
    #                                    , tmpDeepGMLabelMask
    #                                    , slicer.util.getNode("BabyBrain_regMNI2Native_0GenericAffine")
    #                                    , slicer.util.getNode("BabyBrain_regMNI2Native_1Warp")
    #                                    , True)
    #
    #   self.combineLabels(tmpDeepGMLabelMask, outputVolume, outputVolume)
    #   slicer.util.showStatusMessage("Basal ganglia segmentation is finished...")
    #

    # # if not debugMode: # TODO Remover ou nao o debug mode no modulo...???
    # ######################################################################################
    # # Step  - Cleaning temporaty data (Debug mode: Off)
    # ######################################################################################
    # if jumpedAtlasPropagation:
    #   slicer.mrmlScene.RemoveNode(brainAtlasNode)
    #
    # slicer.mrmlScene.RemoveNode(tmpResampledInputNode)
    # slicer.mrmlScene.RemoveNode(brainMaskNode)
    # slicer.mrmlScene.RemoveNode(tmpBrainMask)
    # slicer.mrmlScene.RemoveNode(tmpCerebellumMask)
    # slicer.mrmlScene.RemoveNode(cerebellumOnlyLabelMask)
    # slicer.mrmlScene.RemoveNode(cerebellumMaskNode)
    # slicer.mrmlScene.RemoveNode(tmpBrainstemMask)
    # slicer.mrmlScene.RemoveNode(brainstemMaskNode)
    # slicer.mrmlScene.RemoveNode(venctriculesMaskNode)
    # slicer.mrmlScene.RemoveNode(tmpVentriculesLabelMask)
    # slicer.mrmlScene.RemoveNode(tmpBrainOnlyNode)
    # slicer.mrmlScene.RemoveNode(tmpCerebellumOnlyVolumeNode)
    # slicer.mrmlScene.RemoveNode(tmpBrainstemOnlyVolumeNode)
    # slicer.mrmlScene.RemoveNode(brainstemOnlyLabelMask)
    # slicer.mrmlScene.RemoveNode(brainstemPlusCerebellumLabelMask)
    # slicer.mrmlScene.RemoveNode(brainOnlyLabelMask)
    # slicer.mrmlScene.RemoveNode(hemispheresMaskNode)
    # slicer.mrmlScene.RemoveNode(tmpHemispheresLabelMask)
    # slicer.mrmlScene.RemoveNode(brainOnlyHemispheresLabelMask)
    # slicer.mrmlScene.RemoveNode(fullBrainSegmentationLabelMask)
    #
    # if estimateBasalGanglia:
    #   slicer.mrmlScene.RemoveNode(tmpDeepGMLabelMask)

    logging.info('Processing completed')
    slicer.util.showStatusMessage("Baby Brain Segmentation is finished")

    return True

  #
  # Atlas Propagation
  #
  def atlasPropagation(self, registrationAlgorithm
                       , fixedNode
                       , movingNode
                       , transformLinear
                       , transformElastic
                       , interpolation
                       , sampling
                       , splineGrid
                       , initMethod
                       , numberOfCores
                       , useQuickRegistration):

    if registrationAlgorithm == "BRAINSFit":
      # Applying the first registration level - Linear (Affine)
      regParams = {}
      regParams["fixedVolume"] = fixedNode.GetID()
      regParams["movingVolume"] = movingNode.GetID()
      regParams["outputVolume"] = movingNode.GetID()
      regParams["samplingPercentage"] = sampling
      regParams["linearTransform"] = transformLinear.GetID()
      regParams["initializeTransformMode"] = initMethod
      regParams["useRigid"] = True
      regParams["useAffine"] = True
      regParams["interpolationMode"] = interpolation

      slicer.cli.run(slicer.modules.brainsfit, None, regParams, wait_for_completion=True)

      # Applying the second registration level - Elastic (Spline)
      regParams = {}
      regParams["fixedVolume"] = fixedNode.GetID()
      regParams["movingVolume"] = movingNode.GetID()
      regParams["samplingPercentage"] = sampling
      regParams["bsplineTransform"] = transformElastic.GetID()
      # regParams['initialTransform'] = transformLinear.GetID()
      regParams["initializeTransformMode"] = "Off"
      regParams["splineGridSize"] = splineGrid
      regParams["useBSpline"] = True
      regParams["interpolationMode"] = interpolation

      slicer.cli.run(slicer.modules.brainsfit, None, regParams, wait_for_completion=True)
    else:
      # Create scripts calling. Since the ANTs tools are only provided to Unix systems, the path pattern is fixed.
      home = expanduser("~")
      # Creating temporary folder in home directory
      os.system("mkdir "+ home +"/tmpANTsBabyBrainSegmentation")
      tmpFolder = home + "/tmpANTsBabyBrainSegmentation"

      # Saving the subject image
      slicer.util.saveNode(fixedNode, tmpFolder + '/subject.nii.gz')
      # Saving the brain template
      slicer.util.saveNode(movingNode, tmpFolder + '/template.nii.gz')

      # Use ANTs registration
      if useQuickRegistration:
        os.system("antsRegistrationSyNQuick.sh -d 3 -f " + tmpFolder + "/subject.nii.gz -m " + tmpFolder + "/template.nii.gz -o " + tmpFolder +"/BabyBrain_regMNI2Native_ -n " + str(numberOfCores))
      else:
        os.system("antsRegistrationSyN.sh -d 3 -f " + tmpFolder + "/subject.nii.gz -m " + tmpFolder + "/template.nii.gz -o " + tmpFolder +"/BabyBrain_regMNI2Native_ -n " + str(numberOfCores))

      # Reading registration tranforms
      (read, regTemplate1Warp) = slicer.util.loadTransform(tmpFolder + '/BabyBrain_regMNI2Native_1Warp.nii.gz', True) # TODO Acertar a leitura das transformadas para nao depender somente do nome... do jeito que ta o modulo nao pode rodar duas vezes seguidas...
      regTemplate1Warp.SetName("BabyBrain_regMNI2Native_1Warp") # brain template to native space (SyN)
      (read, regTemplate0GenericAffine) = slicer.util.loadTransform(tmpFolder + '/BabyBrain_regMNI2Native_0GenericAffine.mat', True)
      regTemplate0GenericAffine.SetName("BabyBrain_regMNI2Native_0GenericAffine")# brain template to native space (affine)

      # Removing files from the modules path
      os.system("rm -R " + tmpFolder)

  #
  # Registration Transform Application
  #
  def applyRegistrationTransforms(self, registrationAlgorithm
                                  , inputVolume
                                  , referenceVolume
                                  , outputVolume
                                  , transformLinear
                                  , transformWarp
                                  , isLabel):

    params = {}
    params["inputVolume"] = inputVolume.GetID()
    params["referenceVolume"] = referenceVolume.GetID()
    params["outputVolume"] = outputVolume.GetID()
    params["warpTransform"] = transformLinear.GetID()
    params["inverseTransform"] = False
    if isLabel:
      params["interpolationMode"] = "NearestNeighbor"
      params["pixelType"] = "binary"
    else:
      params["interpolationMode"] = "Linear"
      params["pixelType"] = "float"

    slicer.cli.run(slicer.modules.brainsresample, None, params, wait_for_completion=True)

    params = {}
    params["inputVolume"] = outputVolume.GetID()
    params["referenceVolume"] = inputVolume.GetID()
    params["outputVolume"] = outputVolume.GetID()
    params["warpTransform"] = transformWarp.GetID()
    params["inverseTransform"] = False
    if isLabel:
      params["interpolationMode"] = "NearestNeighbor"
      params["pixelType"] = "binary"
    else:
      params["interpolationMode"] = "Linear"
      params["pixelType"] = "float"

    slicer.cli.run(slicer.modules.brainsresample, None, params, wait_for_completion=True)

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
  # Generic Brain Tissue Segmentation
  #
  def segmentingTissues(self, inputVolume
                           , imageModality
                           , numberOfTissues
                           , outputLabel
                           , oneTissue = False
                           , inputPriorsFile=None
                           , tissueValue = 2
                           , definitions = [1,2,3]):

    params = {}
    params['inputVolume'] = inputVolume.GetID()
    if inputPriorsFile is not "":
        params['inputPriorsFile'] = inputPriorsFile
    params['imageModality'] = imageModality
    params['numberOfTissues'] = numberOfTissues
    params['outputLabel'] = outputLabel.GetID()
    params['oneTissue'] = oneTissue
    params['typeTissue'] = tissueValue
    params['labelsDefinition'] = definitions

    slicer.cli.run(slicer.modules.bayesiantissueclassifier, None, params, wait_for_completion=True)

  #
  # Combine Labels
  #
  def combineLabels(self, firstLabel
                    , secondLabel
                    , outputLabel
                    , firstOverwrites = True):

    params = {}
    params['InputLabelMap_A'] = firstLabel
    params['InputLabelMap_B'] = secondLabel
    params['OutputLabelMap'] = outputLabel
    params['FirstOverwrites'] = firstOverwrites

    slicer.cli.run(slicer.modules.imagelabelcombine, None, params, wait_for_completion=True)


class BabyBrainSegmentationTest(ScriptedLoadableModuleTest):
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
    self.test_BabyBrainSegmentation1()

  def test_BabyBrainSegmentation1(self):
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
    logic = BabyBrainSegmentationLogic()
    self.assertIsNotNone( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')
