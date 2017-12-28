import os
import unittest
import vtk, qt, ctk, slicer
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
    self.parent.title = "Baby Brain Preparation" # TODO make this more human readable by adding spaces
    self.parent.categories = ["Examples"]
    self.parent.dependencies = []
    self.parent.contributors = ["Antonio Carlos da S. Senra Filho (University of Sao Paulo) and Sara"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
This module offers a set of algorithms to brain tissue preparation ... in structural MRI images, namely T1, T2 and PD. This module is the main function to ..., being at moment the naive Bayes classifier available. More details are found in the wikipage: https://www.slicer.org/wiki/Documentation/Nightly/Modules/BabyBrainPreparation
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
    parametersCollapsibleButton.text = "Parameters"
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
    self.inputSelector.noneEnabled = False
    self.inputSelector.showHidden = False
    self.inputSelector.showChildNodeTypes = False
    self.inputSelector.setMRMLScene( slicer.mrmlScene )
    self.inputSelector.setToolTip( "Pick the input to the algorithm. This should be an MRI strutural images with a type listed in the Image Modality option." )
    parametersFormLayout.addRow("Input Volume: ", self.inputSelector)

    #
    # Image Modality
    #
    self.setImageModalityBooleanWidget = ctk.ctkComboBox()
    self.setImageModalityBooleanWidget.addItem("T1")
    self.setImageModalityBooleanWidget.addItem("T2")
    self.setImageModalityBooleanWidget.addItem("PD")
    self.setImageModalityBooleanWidget.setToolTip(
      "MRI strutural image inserted as a input volume.")
    parametersFormLayout.addRow("Image Modality ", self.setImageModalityBooleanWidget)


    #
    # output volume selector TODO Add rename option 
    #
    self.outputSelector = slicer.qMRMLNodeComboBox()
    self.outputSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.outputSelector.selectNodeUponCreation = True
    self.outputSelector.addEnabled = True 
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
    self.setImageResamplingBooleanWidget = ctk.ctkCheckBox()
    self.setImageResamplingBooleanWidget.setChecked(True)
    self.setImageResamplingBooleanWidget.setToolTip(
      "Apply an image space resampling in the input image in order to reduce the total processing time. This is helpful for image with large scale space (usually higher than 256^3 mm).")
    parametersEnhancementLayout.addRow("Apply Image Space Resampling", self.setImageResamplingBooleanWidget)

    #
    # Apply AAD filtering
    #
    self.setApplyAADBooleanWidget = ctk.ctkCheckBox()
    self.setApplyAADBooleanWidget.setChecked(True)
    self.setApplyAADBooleanWidget.setToolTip(
      "Apply the AAD filter on the input data. This is recommended because the image noise level may affect the segmentation performance.")
    parametersEnhancementLayout.addRow("Apply AAD filter", self.setApplyAADBooleanWidget)


    #
    # Apply Bias Correction
    #
    self.setApplyBiasCorrectionBooleanWidget = ctk.ctkCheckBox()
    self.setApplyBiasCorrectionBooleanWidget.setChecked(True)
    self.setApplyBiasCorrectionBooleanWidget.setToolTip(
      "Apply a bias field correction in the input data. This is recommended because the global signal fluctuation provided by magnetic field inhomogeneity may affect the segmentation performance.")
    parametersEnhancementLayout.addRow("Apply Bias Field Correction", self.setApplyBiasCorrectionBooleanWidget)

    #
    # Apply Global Contrast Enhancement Function
    #
    self.setApplyGlobalEnhancementBooleanWidget = ctk.ctkCheckBox()
    self.setApplyGlobalEnhancementBooleanWidget.setChecked(True)
    self.setApplyGlobalEnhancementBooleanWidget.setToolTip(
      "Apply a full image space contrast enhancement in order to increase tissues delineation. This is recommended because the original image may present a poor separation between brain tissues.")
    parametersEnhancementLayout.addRow("Apply Global Contrast Enhancement Function", self.setApplyGlobalEnhancementBooleanWidget)


    #
    # Image Space Resampling Parameters Area
    #
    parametersImageResamplingCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersImageResamplingCollapsibleButton.text = "Image Space Resampling Parameters"
    self.layout.addWidget(parametersImageResamplingCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersImageResamplingLayout = qt.QFormLayout(parametersImageResamplingCollapsibleButton)

    # TODO Colocar uma entrada vector float para o ResampleScalarVolume CLI module -- similar ao N4ITK

    #
    # Interpolation Functions
    #
    self.setInterpolationFunctionBooleanWidget = ctk.ctkComboBox()
    self.setInterpolationFunctionBooleanWidget.addItem("bspline")
    self.setInterpolationFunctionBooleanWidget.addItem("linear")
    self.setInterpolationFunctionBooleanWidget.addItem("nearestNeighbor")
    self.setInterpolationFunctionBooleanWidget.addItem("hamming")
    self.setInterpolationFunctionBooleanWidget.addItem("cosine")
    self.setInterpolationFunctionBooleanWidget.addItem("welch")
    self.setInterpolationFunctionBooleanWidget.addItem("lanczos")
    self.setInterpolationFunctionBooleanWidget.addItem("blackman")
    self.setInterpolationFunctionBooleanWidget.setToolTip(
      "Interpolation functions.")
    parametersImageResamplingLayout.addRow("Interpolation ", self.setInterpolationFunctionBooleanWidget)

    #
    # Noise Attenuation Parameters Area
    #
    parametersNoiseAttenuationCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersNoiseAttenuationCollapsibleButton.text = "Noise Attenuation Parameters"
    self.layout.addWidget(parametersNoiseAttenuationCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersNoiseAttenuationLayout = qt.QFormLayout(parametersNoiseAttenuationCollapsibleButton)

    #
    # Apply Automatic Conductance Function
    #
    self.setApplyAutoConductanceBooleanWidget = ctk.ctkCheckBox()
    self.setApplyAutoConductanceBooleanWidget.setChecked(True)
    self.setApplyAutoConductanceBooleanWidget.setToolTip(
      "Apply an automatic conductance function in order to adjust the conductance parameters using the input image as source. This is recommended because the tissues borders are usually difficult to be empirically infered. If this option is checked, the manual conductance adjustment will be neglected.")
    parametersNoiseAttenuationLayout.addRow("Apply Automatic Conductance Function", self.setApplyAutoConductanceBooleanWidget)

    #
    # Auto Conductance Functions
    #
    self.setConductanceFunctionBooleanWidget = ctk.ctkComboBox()
    self.setConductanceFunctionBooleanWidget.addItem("Canny")
    self.setConductanceFunctionBooleanWidget.addItem("MAD")
    self.setConductanceFunctionBooleanWidget.addItem("Morphological")
    self.setConductanceFunctionBooleanWidget.setToolTip(
      "Automatic conductance functions.")
    parametersNoiseAttenuationLayout.addRow("Conductance Function ", self.setConductanceFunctionBooleanWidget)

    #
    # Apply Conductance Regularization
    #
    self.setApplyConductanceRegularizationBooleanWidget = ctk.ctkCheckBox()
    self.setApplyConductanceRegularizationBooleanWidget.setChecked(True)
    self.setApplyConductanceRegularizationBooleanWidget.setToolTip(
      "Apply an regularization level to the automatic conductance value obtained above. This is recommended only if the automatic conductance function is overestimating the edge cut level, which may return undesired results (tissue border deterioration). This factor should be interpreted as a reduction level for less conservative edge preservation. This option will be applyed only if the automatic conductance choice is used.")
    parametersNoiseAttenuationLayout.addRow("Apply Conductance Regularization", self.setApplyConductanceRegularizationBooleanWidget)

    #
    # Conductance Regularization Value
    #
    self.setFilteringCondutanceWidget = ctk.ctkSliderWidget()
    self.setFilteringCondutanceWidget.maximum = 5
    self.setFilteringCondutanceWidget.minimum = 1
    self.setFilteringCondutanceWidget.value = 2
    self.setFilteringCondutanceWidget.singleStep = 1
    self.setFilteringCondutanceWidget.setToolTip("Conductance regularization.")
    parametersNoiseAttenuationLayout.addRow("Conductance Regularization", self.setFilteringCondutanceWidget)

    #
    # Filtering Parameters: Condutance TODO Adicionar um separador para nao ficar muito junto dos outros paramatros...talvez pular uma linha ou somente colocar um subtitulo "Manual Conductance Adjustment"
    #
    self.setFilteringCondutanceWidget = ctk.ctkSliderWidget()
    self.setFilteringCondutanceWidget.maximum = 50
    self.setFilteringCondutanceWidget.minimum = 0
    self.setFilteringCondutanceWidget.value = 10
    self.setFilteringCondutanceWidget.singleStep = 0.1
    self.setFilteringCondutanceWidget.setToolTip("Conductance parameter.")
    parametersNoiseAttenuationLayout.addRow("Conductance ", self.setFilteringCondutanceWidget)

    #
    # Filtering Parameters: Number of iterations
    #
    self.setFilteringNumberOfIterationWidget = ctk.ctkSliderWidget()
    self.setFilteringNumberOfIterationWidget.maximum = 50
    self.setFilteringNumberOfIterationWidget.minimum = 1
    self.setFilteringNumberOfIterationWidget.value = 25
    self.setFilteringNumberOfIterationWidget.singleStep = 1
    self.setFilteringNumberOfIterationWidget.setToolTip("Number of iterations parameter.")
    parametersNoiseAttenuationLayout.addRow("Number Of Iterations ", self.setFilteringNumberOfIterationWidget)

    #
    # Filtering Parameters: Q value
    #
    self.setFilteringQWidget = ctk.ctkSliderWidget()
    self.setFilteringQWidget.singleStep = 0.01
    self.setFilteringQWidget.minimum = 0.01
    self.setFilteringQWidget.maximum = 1.99
    self.setFilteringQWidget.value = 1.4
    self.setFilteringQWidget.setToolTip("Q value parameter.")
    parametersNoiseAttenuationLayout.addRow("Q Value ", self.setFilteringQWidget)


    #
    # Bias Field Correction Parameters Area
    #
    parametersBiasCorrectionCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersBiasCorrectionCollapsibleButton.text = "Bias Field Correction Parameters"
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
    self.inputMaskSelector.setToolTip( "Pick the input brain mask. This image mask should informs the region of interest where the bias field correction will be applied. If not provided, an Otsu segmentation is performed." )
    parametersBiasCorrectionLayout.addRow("Image Mask: ", self.inputMaskSelector)

    # TODO Colocar um entrada vector float para o N4ITK -- duas entradas para 1 e 2 level

    #
    # Global Contrast Enhancement Parameters Area
    #
    parametersGlobalContrastEnhancementCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersGlobalContrastEnhancementCollapsibleButton.text = "Global Contrast Enhancement Parameters"
    self.layout.addWidget(parametersGlobalContrastEnhancementCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersGlobalContrastEnhancementLayout = qt.QFormLayout(parametersGlobalContrastEnhancementCollapsibleButton)

    #
    # Contrast Modulation Functions TODO Montar metodos CLI para cada uma destas funcoes...fica interessante deixar estes independentes para usar em outros metodos no SLicer
    #
    self.setContrastModulationFunctionBooleanWidget = ctk.ctkComboBox()
    self.setContrastModulationFunctionBooleanWidget.addItem("Logistic")
    self.setContrastModulationFunctionBooleanWidget.addItem("CLAHE") #TODO Implementar CLAHE e CDF image contrast enhancement functions para ITK!!!
    self.setContrastModulationFunctionBooleanWidget.addItem("CDF")
    self.setContrastModulationFunctionBooleanWidget.setToolTip(
      "Contrast modulation functions used to enhance tissue signal contrast in all image scale space. These methods are based on image histogram.")
    parametersGlobalContrastEnhancementLayout.addRow("Contrast Modulation Function ", self.setContrastModulationFunctionBooleanWidget)

    #
    # threshold value
    #
    self.imageThresholdSliderWidget = ctk.ctkSliderWidget()
    self.imageThresholdSliderWidget.singleStep = 0.1
    self.imageThresholdSliderWidget.minimum = -100
    self.imageThresholdSliderWidget.maximum = 100
    self.imageThresholdSliderWidget.value = 0.5
    self.imageThresholdSliderWidget.setToolTip("Set threshold value for computing the output image. Voxels that have intensities lower than this value will set to zero.")
    parametersFormLayout.addRow("Image threshold", self.imageThresholdSliderWidget)

    #
    # check box to trigger taking screen shots for later use in tutorials
    #
    self.enableScreenshotsFlagCheckBox = qt.QCheckBox()
    self.enableScreenshotsFlagCheckBox.checked = 0
    self.enableScreenshotsFlagCheckBox.setToolTip("If checked, take screen shots for tutorials. Use Save Data to write them to disk.")
    parametersFormLayout.addRow("Enable Screenshots", self.enableScreenshotsFlagCheckBox)

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
    enableScreenshotsFlag = self.enableScreenshotsFlagCheckBox.checked
    imageThreshold = self.imageThresholdSliderWidget.value
    logic.run(self.inputSelector.currentNode(), self.outputSelector.currentNode(), imageThreshold, enableScreenshotsFlag)

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

  def takeScreenshot(self,name,description,type=-1):
    # show the message even if not taking a screen shot
    slicer.util.delayDisplay('Take screenshot: '+description+'.\nResult is available in the Annotations module.', 3000)

    lm = slicer.app.layoutManager()
    # switch on the type to get the requested window
    widget = 0
    if type == slicer.qMRMLScreenShotDialog.FullLayout:
      # full layout
      widget = lm.viewport()
    elif type == slicer.qMRMLScreenShotDialog.ThreeD:
      # just the 3D window
      widget = lm.threeDWidget(0).threeDView()
    elif type == slicer.qMRMLScreenShotDialog.Red:
      # red slice window
      widget = lm.sliceWidget("Red")
    elif type == slicer.qMRMLScreenShotDialog.Yellow:
      # yellow slice window
      widget = lm.sliceWidget("Yellow")
    elif type == slicer.qMRMLScreenShotDialog.Green:
      # green slice window
      widget = lm.sliceWidget("Green")
    else:
      # default to using the full window
      widget = slicer.util.mainWindow()
      # reset the type so that the node is set correctly
      type = slicer.qMRMLScreenShotDialog.FullLayout

    # grab and convert to vtk image data
    qimage = ctk.ctkWidgetsUtils.grabWidget(widget)
    imageData = vtk.vtkImageData()
    slicer.qMRMLUtils().qImageToVtkImageData(qimage,imageData)

    annotationLogic = slicer.modules.annotations.logic()
    annotationLogic.CreateSnapShot(name, description, type, 1, imageData)

  def run(self, inputVolume, outputVolume, imageThreshold, enableScreenshots=0):
    """
    Run the actual algorithm
    """

    if not self.isValidInputOutputData(inputVolume, outputVolume):
      slicer.util.errorDisplay('Input volume is the same as output volume. Choose a different output volume.')
      return False

    logging.info('Processing started')

    # Compute the thresholded output volume using the Threshold Scalar Volume CLI module
    cliParams = {'InputVolume': inputVolume.GetID(), 'OutputVolume': outputVolume.GetID(), 'ThresholdValue' : imageThreshold, 'ThresholdType' : 'Above'}
    cliNode = slicer.cli.run(slicer.modules.thresholdscalarvolume, None, cliParams, wait_for_completion=True)

    # Capture screenshot
    if enableScreenshots:
      self.takeScreenshot('BabyBrainPreparationTest-Start','MyScreenshot',-1)

    logging.info('Processing completed')

    return True


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
