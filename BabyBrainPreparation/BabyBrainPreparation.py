import os
import multiprocessing
import platform
import sys
import unittest
import vtk, qt, ctk, slicer
import SimpleITK as sitk
import sitkUtils
from slicer.ScriptedLoadableModule import *
import logging

from os.path import expanduser
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
    # Image Modality
    #
    self.setImageModalityBooleanWidget = ctk.ctkComboBox()
    self.setImageModalityBooleanWidget.addItem("T2")
    self.setImageModalityBooleanWidget.addItem("T1")
    self.setImageModalityBooleanWidget.setToolTip(
      "MRI strutural image inserted as a input volume.")
    parametersFormLayout.addRow("Image Modality ", self.setImageModalityBooleanWidget)

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

    # #
    # # Apply Image Space Resampling
    # #
    # self.setApplyImageResamplingBooleanWidget = ctk.ctkCheckBox()
    # self.setApplyImageResamplingBooleanWidget.setChecked(True)
    # self.setApplyImageResamplingBooleanWidget.setToolTip(
    #   "Apply an image space resampling in the input image in order to reduce the total processing time. This is helpful for "
    #   "image with large scale space (usually higher than 256^3 mm).")
    # parametersEnhancementLayout.addRow("Apply Image Space Resampling", self.setApplyImageResamplingBooleanWidget)

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
    parametersImageResamplingCollapsibleButton.collapsed = False
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
    # Atlas Propagation Parameters Area
    #
    parametersAtlasPropagationCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersAtlasPropagationCollapsibleButton.text = "Atlas Propagation Parameters"
    parametersAtlasPropagationCollapsibleButton.collapsed = True
    self.layout.addWidget(parametersAtlasPropagationCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersAtlasPropagationLayout = qt.QFormLayout(parametersAtlasPropagationCollapsibleButton)

    #
    # Apply Brain Volume Refinement
    #
    self.setApplyBrainVolumeRefinementBooleanWidget = ctk.ctkCheckBox()
    self.setApplyBrainVolumeRefinementBooleanWidget.setChecked(True)
    self.setApplyBrainVolumeRefinementBooleanWidget.setToolTip(
        "Check this if you want to refine the brain volume using the atlas brain mask. This step can be helpful if the "
        "brain extraction procedure left some non-brain tissues in the input image.")
    parametersAtlasPropagationLayout.addRow("Apply Brain Volume Refinement",
                                            self.setApplyBrainVolumeRefinementBooleanWidget)

    #
    # Brain Atlas
    #
    self.setBrainAtlasComboBoxWidget = ctk.ctkComboBox()
    self.setBrainAtlasComboBoxWidget.addItem("NEO2012")  # TODO Ver se usa tambem outro template (2015, http://brain-development.org/brain-atlases/multi-structural-neonatal-brain-atlas/)
    # self.setBrainAtlasComboBoxWidget.addItem("NEO2015") # TODO Novo brain atlas com o mesmo padrao do NEO2012... tem mais detalhes de segmentacao
    self.setBrainAtlasComboBoxWidget.addItem("FET2012")  # TODO Preparar cerebellum e brainstem.
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
    self.setSubjectAgeIntegerWidget.setToolTip(
        "Select the subject's age in weeks. This is only used for neonatal and fetal brain atlases. "
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
    self.setRadioANTs.setToolTip(
        "Use the ANTs SyN diffeomorphic algorithm (recommended). If the ANTs tools are not installed in the user's machine, than this option will be not available.")
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
    modality = self.setImageModalityBooleanWidget.currentText
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
    # useResampling = self.setApplyImageResamplingBooleanWidget.isChecked()
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
              , useAAD
              , useN4ITK
              , useGlobalContrastEnhancement
              , modality
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
          , useAAD
          , useN4ITK
          , useGlobalContrastEnhancement
          , modality
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
    # volumesLogic.CreateScalarVolumeFromVolume(slicer.mrmlScene, outputVolume, inputVolume)
    # clonedInput = volumesLogic.CloneVolume(slicer.mrmlScene, inputVolume, "input_preprocessed")



    # TODO Tentar novamente criar uma barra de progresso : slicer.util.createProgressDialog
    # progressbar = slicer.util.createProgressDialog(parent=slicer.util.mainWindow(),autoClose=False)


    # if useResampling:
    #
    # Step  - Space Resampling
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
      # Step  - Bias field correction - level 1
      #
      self.biasFielCorretion(outputVolume
                             , outputVolume
                             , maskVolume
                             , splineGrid_1
                             , shrinkFactor_1
                             , True
                             , voxelResampling)
      slicer.util.showStatusMessage("Bias field correction (first level) is finished...")

    # Setting the standard paths to extension resources folder
    modulePath = os.path.dirname(slicer.modules.babybrainsegmentation.path)
    if platform.system() is "Windows":
      databasePath = modulePath + "\\Resources\\atlases"
    else:
      databasePath = modulePath + "/Resources/atlases"

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
                                                             "\\templates\\template_" + modality + "_" + str(
        setAge) + ".nii.gz", readingParameters, True)
    else:
      readingParameters = {}
      readingParameters['name'] = "brain_template"
      readingParameters['center'] = True
      readingParameters['show'] = False
      (readSuccess, brainAtlasNode) = slicer.util.loadVolume(databasePath +
                                                             "/" + brainAtlas +
                                                             "/templates/template_" + modality + "_" + str(
        setAge) + ".nii.gz", readingParameters, True)

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
                          , outputVolume
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

    slicer.util.showStatusMessage("Atlas propagation is finished...")

    ######################################################################################
    # Step  - Remove Cerebellum and brainstem from the input volume
    ######################################################################################
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
                                                                 "\\cerebellum\\cerebellum_" + str(setAge) + ".nii.gz",
                                                                 readingParameters,
                                                                 True)
    else:
      readingParameters = {}
      readingParameters['name'] = "cerebellum_template_mask"
      readingParameters['center'] = True
      readingParameters['show'] = False
      readingParameters['labelmap'] = True
      (readSuccess, cerebellumMaskNode) = slicer.util.loadVolume(databasePath +
                                                                 "/" + brainAtlas +
                                                                 "/cerebellum/cerebellum_" + str(setAge) + ".nii.gz",
                                                                 readingParameters,
                                                                 True)

    self.applyRegistrationTransforms(registrationAlgorithm
                                     , cerebellumMaskNode
                                     , outputVolume
                                     , tmpCerebellumMask
                                     , slicer.util.getNode("BabyBrain_regMNI2Native_0GenericAffine")
                                     , slicer.util.getNode("BabyBrain_regMNI2Native_1Warp")
                                     , True)

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
                                                                "\\brainstem\\brainstem_" + str(setAge) + ".nii.gz",
                                                                readingParameters,
                                                                True)
    else:
      readingParameters = {}
      readingParameters['name'] = "brainstem_template_mask"
      readingParameters['center'] = True
      readingParameters['show'] = False
      readingParameters['labelmap'] = True
      (readSuccess, brainstemMaskNode) = slicer.util.loadVolume(databasePath +
                                                                "/" + brainAtlas +
                                                                "/brainstem/brainstem_" + str(setAge) + ".nii.gz",
                                                                readingParameters,
                                                                True)

    self.applyRegistrationTransforms(registrationAlgorithm
                                     , brainstemMaskNode
                                     , outputVolume
                                     , tmpBrainstemMask
                                     , slicer.util.getNode("BabyBrain_regMNI2Native_0GenericAffine")
                                     , slicer.util.getNode("BabyBrain_regMNI2Native_1Warp")
                                     , True)

    # Removing brainstem and cerebellum from the input image
    inputImage = sitkUtils.PullVolumeFromSlicer(outputVolume)
    brainstemMaskImage = sitkUtils.PullVolumeFromSlicer(tmpBrainstemMask)
    brainstemMaskImage = sitk.Cast(brainstemMaskImage, inputImage.GetPixelID())  # This is necessary since the MaskNegated filter requires the same pixel type in both images.
    cerebellumMaskImage = sitkUtils.PullVolumeFromSlicer(tmpCerebellumMask)
    cerebellumMaskImage = sitk.Cast(cerebellumMaskImage, inputImage.GetPixelID())  # This is necessary since the MaskNegated filter requires the same pixel type in both images.
    filter = sitk.MaskNegatedImageFilter()
    output_brainOnly_Image = filter.Execute(inputImage, brainstemMaskImage)
    output_brainOnly_Image = filter.Execute(output_brainOnly_Image, cerebellumMaskImage)

    tmpBrainOnlyNode = slicer.vtkMRMLScalarVolumeNode()
    tmpBrainOnlyNode.SetName("resampled_input_brainOnly")
    slicer.mrmlScene.AddNode(tmpBrainOnlyNode)
    sitkUtils.PushVolumeToSlicer(output_brainOnly_Image, tmpBrainOnlyNode)  # TODO Ver como carregar estas labels sem atualizar a scene... nao fica legal ver essas labels intermediarias durante o processamento

    slicer.util.showStatusMessage("Brainstem and cerebellum removal is finished...")

    # Obtaining only the brainstem and cerebellum volume together
    tmpBrainstemPlusCerebellumLabelOnlyNode = slicer.vtkMRMLLabelMapVolumeNode()
    tmpBrainstemPlusCerebellumLabelOnlyNode.SetName("resampled_input_BrainstemPlusCerebellumOnly_mask")
    slicer.mrmlScene.AddNode(tmpBrainstemPlusCerebellumLabelOnlyNode)
    self.combineLabels(tmpBrainstemMask,tmpCerebellumMask,tmpBrainstemPlusCerebellumLabelOnlyNode)

    inputImage = sitkUtils.PullVolumeFromSlicer(outputVolume)
    brainstemCerebellumMaskImage = sitkUtils.PullVolumeFromSlicer(tmpBrainstemPlusCerebellumLabelOnlyNode)
    # brainstemCerebellumMaskImage = sitk.Cast(brainstemCerebellumMaskImage, inputImage.GetPixelID())  # This is necessary since the MaskNegated filter requires the same pixel type in both images.
    filter = sitk.MaskImageFilter()
    output_brainstemCerebellumOnly_Image = filter.Execute(inputImage, brainstemCerebellumMaskImage)

    tmpBrainstemCerebellumOnlyNode = slicer.vtkMRMLScalarVolumeNode()
    tmpBrainstemCerebellumOnlyNode.SetName("resampled_input_BrainstemPlusCerebellumOnly")
    slicer.mrmlScene.AddNode(tmpBrainstemCerebellumOnlyNode)
    sitkUtils.PushVolumeToSlicer(output_brainstemCerebellumOnly_Image, tmpBrainstemCerebellumOnlyNode)


    if useAAD:
      # Applying noise attenuation on Brainstem+Cerebellum image part. Due to small region, the parameters were adjusted following a empirical analysis.
      #
      # Step  - Noise attenuation
      #
      self.anisotropicAnomalousDiffusionFilter(tmpBrainstemCerebellumOnlyNode
                                               , tmpBrainstemCerebellumOnlyNode
                                               , 10
                                               , 1.2
                                               , timeStep
                                               , conductance
                                               , True
                                               , "Morphological"
                                               , True
                                               , 0.125)
      slicer.util.showStatusMessage("Image noise attenuation (brainstem+cerebellum region) is finished...")

      # Applying noise attenuation on brain image part
      #
      # Step  - Noise attenuation
      #
      self.anisotropicAnomalousDiffusionFilter(tmpBrainOnlyNode
                                               , tmpBrainOnlyNode
                                               , numInt
                                               , qValue
                                               , timeStep
                                               , conductance
                                               , applyAutoConductance
                                               , conductanceFunction
                                               , applyConductanceRegularization
                                               , conductanceRegularization)
      slicer.util.showStatusMessage("Image noise attenuation (brain hemispheres region) is finished...")
      # progressbar.value = 40
      # progressbar.labelText = "Image noise attenuation..."

    if useN4ITK:
      #
      # Step  - Bias field correction - level 2 - Brainstem+Cerebellum
      #
      self.biasFielCorretion(tmpBrainstemCerebellumOnlyNode
                             , tmpBrainstemCerebellumOnlyNode
                             , tmpBrainstemPlusCerebellumLabelOnlyNode
                             , splineGrid_2
                             , shrinkFactor_2
                             , True
                             , voxelResampling)
      slicer.util.showStatusMessage("Bias field correction (second level for brainstem and cerebellum regions) is finished...")

      #
      # Step  - Bias field correction - level 2 - Brain hemispheres
      #
      self.biasFielCorretion(tmpBrainOnlyNode
                             , tmpBrainOnlyNode
                             , maskVolume
                             , splineGrid_2
                             , shrinkFactor_2
                             , True
                             , voxelResampling)
      slicer.util.showStatusMessage(
          "Bias field correction (second level for brain hemispheres regions) is finished...")

    # Adding both parts to the same volume
    self.addVolumesTogether(tmpBrainOnlyNode, tmpBrainstemCerebellumOnlyNode, outputVolume)

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

    ######################################################################################
    # Step  - Cleaning temporaty data (Debug mode: Off)
    ######################################################################################
    slicer.mrmlScene.RemoveNode(brainAtlasNode)
    slicer.mrmlScene.RemoveNode(tmpCerebellumMask)
    slicer.mrmlScene.RemoveNode(cerebellumMaskNode)
    slicer.mrmlScene.RemoveNode(tmpBrainstemMask)
    slicer.mrmlScene.RemoveNode(brainstemMaskNode)
    slicer.mrmlScene.RemoveNode(tmpBrainOnlyNode)


    home = expanduser("~")
    tmpFolder = home + "/tmpANTsBabyBrainSegmentation"
    os.system("rm -R " + tmpFolder)

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
      os.system("mkdir " + home + "/tmpANTsBabyBrainSegmentation")
      tmpFolder = home + "/tmpANTsBabyBrainSegmentation"

      # Saving the subject image
      slicer.util.saveNode(fixedNode, tmpFolder + '/subject.nii.gz')
      # Saving the brain template
      slicer.util.saveNode(movingNode, tmpFolder + '/template.nii.gz')

      # Use ANTs registration
      if useQuickRegistration:
        os.system(
          "antsRegistrationSyNQuick.sh -d 3 -f " + tmpFolder + "/subject.nii.gz -m " + tmpFolder + "/template.nii.gz -o " + tmpFolder + "/BabyBrain_regMNI2Native_ -n " + str(
            numberOfCores))
      else:
        os.system(
          "antsRegistrationSyN.sh -d 3 -f " + tmpFolder + "/subject.nii.gz -m " + tmpFolder + "/template.nii.gz -o " + tmpFolder + "/BabyBrain_regMNI2Native_ -n " + str(
            numberOfCores))

      # Reading registration tranforms
      (read, regTemplate1Warp) = slicer.util.loadTransform(tmpFolder + '/BabyBrain_regMNI2Native_1Warp.nii.gz', True)
      regTemplate1Warp.SetName("BabyBrain_regMNI2Native_1Warp")  # brain template to native space (SyN)
      (read, regTemplate0GenericAffine) = slicer.util.loadTransform(tmpFolder + '/BabyBrain_regMNI2Native_0GenericAffine.mat', True)
      regTemplate0GenericAffine.SetName("BabyBrain_regMNI2Native_0GenericAffine")  # brain template to native space (affine)

      # Removing files from the modules path
      # os.system("rm -R " + tmpFolder)

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
      params["interpolationMode"] = "BSpline"
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
      params["interpolationMode"] = "BSpline"
      params["pixelType"] = "float"

    slicer.cli.run(slicer.modules.brainsresample, None, params, wait_for_completion=True)

  #
  # Combine Labels
  #
  def combineLabels(self, firstLabel
                    , secondLabel
                    , outputLabel
                    , firstOverwrites=True):

    params = {}
    params['InputLabelMap_A'] = firstLabel
    params['InputLabelMap_B'] = secondLabel
    params['OutputLabelMap'] = outputLabel
    params['FirstOverwrites'] = firstOverwrites

    slicer.cli.run(slicer.modules.imagelabelcombine, None, params, wait_for_completion=True)

  #
  # Add Volumes
  #
  def addVolumesTogether(self, volume1Node, volume2Node, outputNode, order=0):
    params = {}
    params['inputVolume1'] = volume1Node
    params['inputVolume2'] = volume2Node
    params['outputVolume'] = outputNode
    params['order'] = order

    slicer.cli.run(slicer.modules.addscalarvolumes, None, params, wait_for_completion=True)

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
    params["useConductanceRegularization"] = useConductanceRegularization
    params["kappaRegFactor"] = conductanceRegularization
    params["useAutoConductance"] = useAutoConductance
    params["optFunction"] = conductanceFunction
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
