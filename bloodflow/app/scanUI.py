''' Scan UI for Bloodflow / Stroke Detection '''

# All dependencies are from the conda 'base' environment.

import datetime
import json
import pathlib
import platform
import socket
import subprocess
import sys
import traceback
import time
windowPlatform = False
if platform.win32_ver()[0] == '10':
  import winreg
  import winsound
  windowPlatform = True

from ctypes import CDLL, c_bool, c_char_p, c_double, c_float, c_int, POINTER
if platform.win32_ver()[0] == '10':
  from ctypes import WINFUNCTYPE as FUNCTYPE
else:
  from ctypes import CFUNCTYPE as FUNCTYPE
import random
from PyQt5.QtWidgets import QMessageBox
from PyQt5.QtGui import QPixmap
#from PyQt5.QtCore import pyqtSlot

# Load the DLL and set up expected arg and return types.
# TODO(jfs): Factor this.
scanDLLPath = '../../bazel-bin/bloodflow/plugin/scan.dll'
scanDLL = None

# use scanDLL with added export functions
extended_scanDLL = True
# number of dta samples to collect and plot
plot_array_size = 200 # default 200, 5040 for 3 minutes, update in Cameramanager.h if >1000
camFrameRate = 28 # default 28 BH

try:
  # Import the DLL, and set non-default return types. Functions are those exported by the DLL.
  # (per https://docs.python.org/3/library/ctypes.html)
  scanDLL = CDLL(scanDLLPath)
  scanDLL.init.restype = c_bool
  scanDLL.start.restype = c_bool
  scanDLL.stop.restype = c_bool
  scanDLL.triggerDataCollection.restype = c_bool
  scanDLL.rcam_model.restype = c_char_p
  scanDLL.waitDataCollection.restype = c_bool
  scanDLL.writeCSV.restype = c_bool
  scanDLL.close.restype = c_bool
  scanDLL.checkData.restype = c_bool
  scanDLL.dataMean.restype = c_float
  scanDLL.dataContrast.restype = c_float
  scanDLL.meanMean.restype = c_double
  scanDLL.darkMean.restype = c_double
  scanDLL.lightMean.restype = c_double
  scanDLL.framesAccumulated.restype = c_bool
  scanDLL.writeAverageFrames.restype = c_bool
  scanDLL.anyBadFrames.restype = c_bool
  if extended_scanDLL:
    scanDLL.dataSTD.restype = c_float
    scanDLL.sampleContrast.restype = c_float
    scanDLL.sampleMean.restype = c_float
    scanDLL.sampleSTD.restype = c_float
    getArrContrast = scanDLL.arrContrast
    getArrContrast.restype = POINTER(c_float * plot_array_size)
    getArrMean = scanDLL.arrMean
    getArrMean.restype = POINTER(c_float * plot_array_size)
    getArrSTD = scanDLL.arrSTD
    getArrSTD.restype = POINTER(c_float * plot_array_size)
    #scanDll.arrSTD.restype = POINTER(c_float * 300)
  print('Imported scan DLL from', scanDLLPath)  # No log widget at this point.
except OSError as err:
  print("OS error loading %s: %s" % (scanDLLPath, err))
  # run anyway, for testing (TODO(jfs): only if not on Windows)
except: # pylint: disable=bare-except
  exc_type0, exc_value0, exc_traceback0 = sys.exc_info()
  traceback.print_tb(exc_traceback0, limit=5, file=sys.stdout)

# This must be imported after CDLL, because numpy mangles path searching on Windows
# See:
#   https://docs.conda.io/projects/conda/en/latest/user-guide/troubleshooting.html#numpy-mkl-library-load-failed
#   https://github.com/ContinuumIO/anaconda-issues/issues/10628
#   https://github.com/numpy/numpy/issues/11431
# There are workarounds, as yet untried.

import numpy as np

import matplotlib# pylint: disable=unused-import
from matplotlib.backends.qt_compat import QtCore, QtGui, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import speckle_fcns
import aligner
if windowPlatform:
  import ctlThorITC
#
# Composite widget classes (like ipywidgets)
#

class BoundedFloatText:
  ''' Text field with a label and float semantics '''
  def __init__(self, description, value=0):
    self.label = QtWidgets.QLabel(description)
    self.lineEdit = QtWidgets.QLineEdit(str(value))
    self.layout = QtWidgets.QHBoxLayout()
    self.layout.addWidget(self.label)
    self.layout.addWidget(self.lineEdit)

  def value(self):
    ''' Return the value in the lineEdit widget. '''
    return float(self.lineEdit.text())

class Textarea:
  ''' Multiline text field with a label '''
  def __init__(self, description, value='', lines=1):
    self.label = QtWidgets.QLabel(description)
    if lines > 1:
      self.textEdit = QtWidgets.QPlainTextEdit(value)
    else:
      self.textEdit = QtWidgets.QLineEdit(value)
    self.layout = QtWidgets.QVBoxLayout()
    self.layout.addWidget(self.label)
    self.layout.addWidget(self.textEdit)
    self.layout.addStretch()

  def value(self):
    ''' Return the value in the textEdit widget. '''
    return self.textEdit.toPlainText()

def Button(description=None, callback=None):
  ''' Button with a title and a callback '''
  button = QtWidgets.QPushButton(description)
  button.clicked.connect(callback)
  button.setStyleSheet("background-color: lightgrey;")
  bf = button.font()
  #bf.setBold(True)
  #fps = button.font().pointSize()
  bf.setPointSize(bf.pointSize()*2)
  button.setFont(bf)
  return button

def DataModeButton(description=None, callback=None):
  ''' Button with a title and a callback '''
  button = QtWidgets.QPushButton(description)
  button.clicked.connect(callback)
  button.setStyleSheet("background-color: lightgrey;")
  bf = button.font()
  #bf.setBold(True)
  #fps = button.font().pointSize()
  bf.setPointSize(bf.pointSize()*2)
  button.setFont(bf)
  return button

def makeExitButton(description, callback):
  ''' Button with a title and a callback '''
  button = QtWidgets.QPushButton(description)
  button.clicked.connect(callback)
  button.setStyleSheet("background-color: lightgrey;")
  bf = button.font()
  bf.setPointSize(bf.pointSize()*2)
  button.setFont(bf)
  return button

# plotting data
class PCanvas(FigureCanvasQTAgg):
  ''' Plot canvas holder '''
  def __init__(self, label, ndata, parent=None, width=5, height=4, dpi=100):# pylint: disable=unused-argument
    ''' ... '''
    self.ndata = ndata
    self.xdata = [x / camFrameRate for x in list(range(ndata))] # should be linked to scanDict['cameraParameters']['frameAcquisitionRate_Hz']
    self.ind = 0
    self.nPlotTimes = 0
    fig = Figure(figsize=(width, height), dpi=dpi)
    self.axes = fig.add_subplot(111)
    self.label = label
    # self.axes.legend([label])
    self.axes.set_xlim(0,self.xdata[-1])
    self.axes.set_ylim(0,1)
    #styles = {'color':'b', 'font-size':'12px'}
    #self.axes.set_ylabel(str(self.label))
    self.ydataMean = [0.0 for i in range(ndata)]
    self.ydataContrast = [0.0 for i in range(ndata)]
    self.ydataSTD = [0.0 for i in range(ndata)]
    self.ydataMean_min = 0
    self.ydataMean_max = 0
    self.ydataContrast_min = 0
    self.ydataContrast_max = 0
    self.ydataSTD_min = 0
    self.ydataSTD_max = 0
    super(PCanvas, self).__init__(fig)

  def clear(self):
    '''clean data'''
    self.ind = 0
    self.nPlotTimes = 0
    self.ydataMean = [0.0 for i in range(self.ndata)]
    self.ydataContrast = [0.0 for i in range(self.ndata)]
    self.ydataSTD = [0.0 for i in range(self.ndata)]

  def dataMean_range(self):
    '''...'''
    if self.nPlotTimes <= 1:
      self.ydataMean_min = 0
      n0 = [i for i in self.ydataMean if i !=0]
      if len(n0) > 0:
        self.ydataMean_min = min(n0)
      self.ydataMean_max = max(self.ydataMean)
      if self.ydataMean_max < 0.1:
        self.ydataMean_max = 120
      if self.ydataMean_max - self.ydataMean_min < 10:
        self.ydataMean_min = self.ydataMean_max * 0.05
  def dataContrast_range(self):
    '''...'''
    if self.nPlotTimes <= 1:
      self.ydataContrast_min = 0
      n0 = [i for i in self.ydataContrast if i != 0]
      if len(n0) > 0:
        self.ydataContrast_min = min(n0)
      self.ydataContrast_max = max(self.ydataContrast)
      if self.ydataContrast_max < 0.1:
        self.ydataContrast_max = 5
      if self.ydataContrast_max - self.ydataContrast_min < 2:
        self.ydataContrast_min = self.ydataContrast_max * 0.2
  def dataSTD_range(self):
    '''...'''
    if self.nPlotTimes <= 1:
      self.ydataSTD_min = 0
      n0 = [i for i in self.ydataSTD if i != 0]
      if len(n0) > 0:
        self.ydataSTD_min = min(n0)
      self.ydataSTD_max = max(self.ydataSTD)
      if self.ydataSTD_max < 0.1:
        self.ydataSTD_max = 5
      if self.ydataSTD_max - self.ydataSTD_min < 2:
        self.ydataSTD_min = self.ydataSTD_max * 0.2
  def plotData(self, dataMode, nskip):
    '''...'''
    self.ind = self.ind + 1
    self.ind = self.ind % self.ndata
    if nskip == 0 or self.ind % nskip == 0:
      self.axes.cla()
      if dataMode == "Contrast":
        # self.axes.set_ylim([self.ydataContrast_min, self.ydataContrast_max*1.1])
        dataHolder_noZeros = [i for i in self.ydataContrast if i != 0]
        if len(dataHolder_noZeros) == 0:
          dataHolder_noZeros = [0, 0]
        self.axes.plot(self.xdata, [x or dataHolder_noZeros[-1] for x in self.ydataContrast], 'b')
      else:
        if dataMode == "Mean":
          # self.axes.set_ylim([self.ydataMean_min, self.ydataMean_max*1.1])
          dataHolder_noZeros = [i for i in self.ydataMean if i != 0]
          if len(dataHolder_noZeros) == 0:
            dataHolder_noZeros = [0, 0]
          self.axes.plot(self.xdata, [x or dataHolder_noZeros[-1] for x in self.ydataMean], 'b')
        else:
          if dataMode == "STD":
            # self.axes.set_ylim([self.ydataSTD_min, self.ydataSTD_max*1.1])
            dataHolder_noZeros = [i for i in self.ydataSTD if i != 0]
            if len(dataHolder_noZeros) == 0:
              dataHolder_noZeros = [0, 0]
            self.axes.plot(self.xdata, [x or dataHolder_noZeros[-1] for x in self.ydataSTD], 'b')
      self.draw()
      QtCore.QCoreApplication.instance().processEvents()
      self.nPlotTimes = self.nPlotTimes + 1
#
# The scanner
#

class Scanner(QtWidgets.QMainWindow):
  ''' Scanner class '''
  def __init__(self, scanDLL_):
    ''' Create a Scanner
    args:
      scanDLL_: camera,etc DLL; injected for testability
    '''
    super().__init__()

    # Check if Thorlab laser drivers are present, if so, connect and turns on TECs
    self.thorLaserAvailable = False # Used to keep non-Thorlab laser systems from attempting control of ThorITCs
    self.thorLaserStatus = False # Prevents control of ThorITCs if any ThorITC-related errors occur
    if windowPlatform:
      self.instrSeed, self.instrTA, self.resMngr, self.thorLaserAvailable = ctlThorITC.openLaserConn()
    if self.thorLaserAvailable:
      self.thorLaserStatus = ctlThorITC.turnOnTec(self.instrSeed,self.instrTA)
    if self.thorLaserAvailable and not self.thorLaserStatus:
      self.log('ERROR: Thorlabs TECs could not be turned on, connection disabled.')
      print('\a')
    # self.thorLaserAvailable = False # Uncomment to override computer controlling of Thorlabs lasers

    # Set up a logging callback for the DLL.
    logFuncType = FUNCTYPE(c_int, c_char_p)
    def logCallback(b):
      self.log(str(b, 'utf-8'))
      return 0
    self.logFunc = logFuncType(logCallback)  # for extra credit: lambda s: self.log(s)

    self.scanDLL = scanDLL_  # DLL with camera, etc C++ code

    # Ok to run (for testing) w/o an octopus?
    self.noOctopusOk_ = False

    # Capture this many images (i.e., image pairs) for each location.
    self.numImages_ = int(plot_array_size / 2)

    # Wait for this long for the trigger to be pressed.
    self.timeout_s_ = 600

    # System is hardware triggered? Assume true, if there's an octopus.
    self.hwTriggered_ = self.scanDLL and self.scanDLL.fx3_NumDevices(0x4F12) == 1

    self.filename_ = 'calibrate'  # ToDo(jfs): UI widget?
    self.saveTestFrames_ = False  # Save dark/light test frames during scan?

    # Start the log file in the machine dir, as we don't have a scanData dir yet.
    # ToDo(jfs): This file grows without bound.
    outputDir = outputDirectory()
    pathlib.Path(outputDir).mkdir(parents=True, exist_ok=True)
    self.logFile = open(outputDir + '/log.txt', 'a')

    # Graphics for capture 12 locations
    # self.locations_ = [
    #   [ 'rsrc/head_L.jpg', 'rsrc/head_R.jpg' ],  # quiescent
    #   [ 'rsrc/head_L1.jpg', 'rsrc/head_R0.jpg' ],  # location 1, etc.
    #   [ 'rsrc/head_L1c.jpg', 'rsrc/head_R1.jpg' ],
    #   [ 'rsrc/head_L2.jpg', 'rsrc/head_R1c.jpg' ],
    #   [ 'rsrc/head_L2c.jpg', 'rsrc/head_R2.jpg' ],
    #   [ 'rsrc/head_L3.jpg', 'rsrc/head_R2c.jpg' ],
    #   [ 'rsrc/head_L3c.jpg', 'rsrc/head_R3.jpg' ],
    #   [ 'rsrc/head_L4.jpg', 'rsrc/head_R3c.jpg' ],
    #   [ 'rsrc/head_L4c.jpg', 'rsrc/head_R4.jpg' ],
    #   [ 'rsrc/head_L5.jpg', 'rsrc/head_R4c.jpg' ],
    #   [ 'rsrc/head_L5c.jpg', 'rsrc/head_R5.jpg' ],
    #   [ 'rsrc/head_L6.jpg', 'rsrc/head_R5c.jpg' ],
    #   [ 'rsrc/head_L6c.jpg', 'rsrc/head_R6.jpg' ],
    #   [ 'rsrc/head_L6c.jpg', 'rsrc/head_R6c.jpg' ]  # done
    # ]
    # Graphics for capture 30 locations
    self.locations_ = [
      [ 'rsrc/00L_Start_Empty.jpg', 'rsrc/00R_Start_Empty.jpg' ], # empty
      [ 'rsrc/11L_Forehead_Vert.jpg', 'rsrc/00R_Start_Empty.jpg' ], # location L1
      [ 'rsrc/12L_Forehead_Vert.jpg', 'rsrc/00R_Start_Empty.jpg' ],
      [ 'rsrc/13L_Forehead_Vert.jpg', 'rsrc/00R_Start_Empty.jpg' ],
      [ 'rsrc/21L_Forehead_Hori.jpg', 'rsrc/00R_Start_Empty.jpg' ],
      [ 'rsrc/22L_Forehead_Hori.jpg', 'rsrc/00R_Start_Empty.jpg' ],
      [ 'rsrc/23L_Forehead_Hori.jpg', 'rsrc/00R_Start_Empty.jpg' ],
      [ 'rsrc/31L_SilvFiss_Angle.jpg', 'rsrc/00R_Start_Empty.jpg' ],
      [ 'rsrc/32L_SilvFiss_Angle.jpg', 'rsrc/00R_Start_Empty.jpg' ],
      [ 'rsrc/33L_SilvFiss_Angle.jpg', 'rsrc/00R_Start_Empty.jpg' ],
      [ 'rsrc/41L_Temple_Vert.jpg', 'rsrc/00R_Start_Empty.jpg' ],
      [ 'rsrc/42L_Temple_Vert.jpg', 'rsrc/00R_Start_Empty.jpg' ],
      [ 'rsrc/43L_Temple_Vert.jpg', 'rsrc/00R_Start_Empty.jpg' ],
      [ 'rsrc/51L_Temple_Hori.jpg', 'rsrc/00R_Start_Empty.jpg' ],
      [ 'rsrc/52L_Temple_Hori.jpg', 'rsrc/00R_Start_Empty.jpg' ],
      [ 'rsrc/53L_Temple_Hori.jpg', 'rsrc/00R_Start_Empty.jpg' ], # location L15
      [ 'rsrc/00L_Start_Empty.jpg', 'rsrc/11R_Forehead_Vert.jpg' ], # location R1
      [ 'rsrc/00L_Start_Empty.jpg', 'rsrc/12R_Forehead_Vert.jpg' ],
      [ 'rsrc/00L_Start_Empty.jpg', 'rsrc/13R_Forehead_Vert.jpg' ],
      [ 'rsrc/00L_Start_Empty.jpg', 'rsrc/21R_Forehead_Hori.jpg' ],
      [ 'rsrc/00L_Start_Empty.jpg', 'rsrc/22R_Forehead_Hori.jpg' ],
      [ 'rsrc/00L_Start_Empty.jpg', 'rsrc/23R_Forehead_Hori.jpg' ],
      [ 'rsrc/00L_Start_Empty.jpg', 'rsrc/31R_SilvFiss_Angle.jpg' ],
      [ 'rsrc/00L_Start_Empty.jpg', 'rsrc/32R_SilvFiss_Angle.jpg' ],
      [ 'rsrc/00L_Start_Empty.jpg', 'rsrc/33R_SilvFiss_Angle.jpg' ],
      [ 'rsrc/00L_Start_Empty.jpg', 'rsrc/41R_Temple_Vert.jpg' ],
      [ 'rsrc/00L_Start_Empty.jpg', 'rsrc/42R_Temple_Vert.jpg' ],
      [ 'rsrc/00L_Start_Empty.jpg', 'rsrc/43R_Temple_Vert.jpg' ],
      [ 'rsrc/00L_Start_Empty.jpg', 'rsrc/51R_Temple_Hori.jpg' ],
      [ 'rsrc/00L_Start_Empty.jpg', 'rsrc/52R_Temple_Hori.jpg' ],
      [ 'rsrc/00L_Start_Empty.jpg', 'rsrc/53R_Temple_Hori.jpg' ], # location R15
      [ 'rsrc/00L_Start_Empty.jpg', 'rsrc/00R_Start_Empty.jpg' ]  # done
    ]

    # Load cameraInfo for source-to-detector distance (and anything else we add).
    with open('cameras.json', 'r') as cameraInfo:
      self._cameras = json.load(cameraInfo)

    # Create the outer layer of the UI. Everything else attches to this.
    self._main = QtWidgets.QWidget()
    self.setCentralWidget(self._main)
    # maind window layout will also include bloodflow data view window
    main_layout = QtWidgets.QHBoxLayout(self._main)
    layout = QtWidgets.QVBoxLayout()

    # logo
    logo = QtWidgets.QLabel()
    logo.setPixmap(QtGui.QPixmap('rsrc/Opernwater_logo_gradient_100.png'))
    layout.addWidget(logo)

    # scale factor for pixmaps relative to a base screen resolution
    screen_sz = QtCore.QCoreApplication.instance().primaryScreen().size()
    base_screen_h = 1080
    self.pixmap_scale = base_screen_h / screen_sz.height()

    # scan location guide images
    hbox = QtWidgets.QHBoxLayout()
    self.picRight = QtWidgets.QLabel()

#    self.picRight.setPixmap(QtGui.QPixmap(self.locations_[0][1]))
    img = QtGui.QPixmap(self.locations_[0][1])
    img_height = img.height()
    self.picRight.setPixmap(img.scaledToHeight(img_height / self.pixmap_scale))
    #self.picRight.setPixmap(img)
    hbox.addWidget(self.picRight)
    self.picLeft = QtWidgets.QLabel()

    img = QtGui.QPixmap(self.locations_[0][0])
    img_height = img.height()
    self.picLeft.setPixmap(img.scaledToHeight(img_height / self.pixmap_scale))
    #self.picLeft.setPixmap(img)
    hbox.addWidget(self.picLeft)
    layout.addLayout(hbox)

    # subject ID and notes field for whatever the diagnostician wants to call out
    vbox = QtWidgets.QVBoxLayout()
    self.subjectID = Textarea('Subject ID', lines=1)
    self.subjectID.label.setMaximumWidth(75)
    self.subjectID.textEdit.setMaximumWidth(200)
    vbox.addLayout(self.subjectID.layout)
    self.experimentNotes = Textarea('Notes (operator, position, repeat #)', lines=3)
    vbox.addLayout(self.experimentNotes.layout)
    layout.addLayout(vbox)

    #container of buttons, except "Exit" button
    self.buttonLock = False
    self.buttons = []
    # scan control widgets
    hbox = QtWidgets.QHBoxLayout()
    if len(sys.argv) >= 2 and sys.argv[1] == '--align':
      self.alignButton = Button('Align', self.onAlignClick)
      hbox.addWidget(self.alignButton)
      self.buttons.append(self.alignButton)

    self.calibrateButton = Button('Calibrate', self.onCalibrateClick)
    self.buttons.append(self.calibrateButton)
    self.scanButton = Button('Scan', self.onScanClick)
    self.buttons.append(self.scanButton)
    self.backupButton = Button('Backup', self.onBackupClick)
    self.buttons.append(self.backupButton)
    self.exitButton = makeExitButton('Exit', self.onExitButtonClick)

    hbox.addWidget(self.calibrateButton)
    hbox.addWidget(self.scanButton)
    hbox.addWidget(self.backupButton)
    hbox.addWidget(self.exitButton)
    layout.addLayout(hbox)

    # log area: all output should go here!
    self.logText = QtWidgets.QPlainTextEdit()
    layout.addWidget(self.logText)

    #main_layout.addLayout(layout)
    self.camerasSN = []

    # bloodflow data view layout
    data_layout = QtWidgets.QVBoxLayout()
    data_layout.setAlignment(QtCore.Qt.AlignTop)
    # view mode buttons
    self.button_contrast = DataModeButton("Contrast", self.onModeContrastClick)
    self.button_mean = DataModeButton("Mean", self.onModeMeanClick)
    self.button_std = DataModeButton("Standard Deviation", self.onModeSTDClick)
    mode_buttons_layout = QtWidgets.QHBoxLayout()
    mode_buttons_layout.addWidget(self.button_contrast)
    mode_buttons_layout.addWidget(self.button_mean)
    mode_buttons_layout.addWidget(self.button_std)
    data_layout.addLayout(mode_buttons_layout)
    # plot area
    self.plot_layout = QtWidgets.QVBoxLayout()
    self.ndata_plot = plot_array_size
    self.plots = []
    self.calibrationFolder = ''
    self.createDataPlotForCameras()

    data_layout.addLayout(self.plot_layout)

    # compose main layout
    main_layout.addLayout(data_layout)
    main_layout.addLayout(layout)

    #default data plot mode is "Contrast"
    self.dataMode = "Contrast"

    self.nplot_updates = 0
    self.nplot_updates_period = 3 #!!!
    self.onModeContrastClick()

    self.clearPlots()

    # Log the software version.
    gitInfo = self.getGitInfo()
    self.log('Version: %s' % gitInfo['tagVersion'])

    self.alignerObj = None

    self.state = "idle"

    #!!! set debug option to False before committing
    self.enable_setting_debug_options = False
    if not self.enable_setting_debug_options:
      # debug options disabled
      self.use_software_trigger_for_debug = False
      self.skip_calibration_for_debug = False
      self.alignWindowModeless = False
      self.testData = False #!!! set this True to watch random data plots in "scan" and "align" modes
      self.testAcceptBadFrames = False
    else:
      #debug options enabled
      self.use_software_trigger_for_debug = True
      self.skip_calibration_for_debug = True
      self.alignWindowModeless = False
      self.testData = False #!!! set this True to watch random data plots in "scan" and "align" modes
      self.testAcceptBadFrames = True

    if self.use_software_trigger_for_debug:
      self.log("!!!!!!!!!!!!!!!!!!!!! Debug option is set: use_software_trigger_for_debug")
    if self.skip_calibration_for_debug:
      self.log("!!!!!!!!!!!!!!!!!!!!! Debug option is set: skip_calibration_for_debug")
    if self.alignWindowModeless:
      self.log("!!!!!!!!!!!!!!!!!!!!! Debug option is set: alignWindowModeless")
    if self.testData:
      self.log("!!!!!!!!!!!!!!!!!!!!! Debug option is set: test data plotting")
    if self.testAcceptBadFrames:
      self.log("!!!!!!!!!!!!!!!!!!!!! Debug option is set: testAcceptBadFrames")

    # Are the cameras calibrated? Requires running aligner first.
    self.calibrated_ = self.skip_calibration_for_debug

  def removeDataPlotsFromLayout(self):
    ''' ... '''
    for i in reversed(range(self.plot_layout.count())):
      item = self.plot_layout.itemAt(i)
      if item:
        w = item.widget()
        if w:
          w.setParent(None)
          w.deleteLater()
    self.plots = []
    return True

  def startCameras(self):
    ''' Start the cameras. '''
    # Create and write scan metadata.
    #filename = 'scan_check'
    scanDict = self.createScanDict()
    scanDataDir = scanDict['fileParameters']['localScanDataDir']
    metadataFileName = scanDataDir + '/scan_metadata.json'
    self.writeScanDict(scanDict, metadataFileName)
    # Point the log to the data output directory for the duration of the scan.
    # Log the scan start time.
    self.log('Scan check: %s' % time.asctime())
    # Initialize the scanner with metadata (creates CameraManager, OctopusManager, etc).
    cFilename = c_char_p(bytes(metadataFileName, 'utf-8'))
    if not self.scanDLL.init(cFilename, self.logFunc):  # This is a re-init, given aligner has run.
      self.log("ERROR: Scanner init failed.")
      return False
    # if not self.scanDLL.start(c_bool(False)): # BH removing because activating laser output at start of UI
    #   self.log("ERROR: Can't start cameras.")
    #   return False
    #if not self.scanDLL.triggerDataCollection():
    #  self.log('ERROR: Data collection trigger failed.')
    #  return False
    return True

  def stopCameras(self):
    ''' Stop the cameras. '''
    if not self.scanDLL.stop():
      self.log('ERROR: ScanDLL stop() failed.')
      return False
    return True

  def setupDataPlotCameras(self):
    ''' Find the attached cameras and set them up. '''
    # Clear any existing image widgets.
    #while self._cameras_layout.count():
    #  w = self._cameras_layout.itemAt(0).widget()
    #  w.setText('')
    #  self._cameras_layout.removeWidget(w)

    # Create a list of cameras
    self.camerasSN = []
    if not self.startCameras():
      self.log("Cameras check failed.")
      return False
    numCameras = self.scanDLL.rcam_NumCameras()
    self.log('Found %d camera(s)' % numCameras)
    for nth in range(numCameras):
      serialNumber = self.scanDLL.rcam_serialNumber(c_int(nth))
      if str(serialNumber) not in self._cameras:
        self.log('"ERROR: Camera %d not in camera.json file; please fix and restart scanner.' % serialNumber )
        self.lock_buttons(True)
        return False
      if serialNumber <= 0:
        self.log("ERROR: Can't open camera index %d" % nth)
      else:
        self.camerasSN.append(serialNumber)
    return True

  def createDataPlotForCameras(self):
    ''' ... '''
    ret = False
    if self.setupDataPlotCameras():
      self.removeDataPlotsFromLayout()
      for cameraSN in self.camerasSN:
        self.plots.append(self.createPlot(cameraSN, self.ndata_plot))
      for p in self.plots:
        self.plot_layout.addWidget(QtWidgets.QLabel(str(p.label)))
        p.setMinimumWidth(600)
        self.plot_layout.addWidget(p)
      ret = True
    return ret

  def clearPlots(self):
    ''' ... '''
    for p in self.plots:
      p.clear()
      p.draw()
    QtCore.QCoreApplication.instance().processEvents()

  def createPlot(self, label, ndata):
    ''' create plot with some default parameters '''
    sc = PCanvas(label, ndata, width=20, height=4, dpi=100)
    return sc

  def update_plot(self):
    ''' Update plot (on Qt GUI timer callback) '''
    try:
      # remove the first, append a new one.
      for p in self.plots:
        cSerialNumber = c_int(p.label)
        if self.testData:
          p.ydataContrast = p.ydataContrast[1:] + [random.randint(0, 10)]
          p.ydataMean = p.ydataMean[1:] + [random.randint(0, 10)]
          p.ydataSTD = p.ydataSTD[1:] + [random.randint(0, 10)]
        else:
          if self.state == "scan" or self.state == "align":
            #p.ydataContrast = getArrContrast(cSerialNumber)[:300]
            i = 0
            for v in getArrContrast(cSerialNumber).contents:
              p.ydataContrast[i] = v
              i = i + 1
            p.dataContrast_range()
            #p.ydataMean = getArrMean(cSerialNumber)[:300]
            i = 0
            for v in getArrMean(cSerialNumber).contents:
              p.ydataMean[i] = v
              i = i + 1
            p.dataMean_range()
            #p.ydataSTD = getArrSTD(cSerialNumber)[:300]
            i = 0
            for v in getArrSTD(cSerialNumber).contents:
              p.ydataSTD[i] = v
              i = i + 1
            p.dataSTD_range()
            #p.ydataContrast[p.ind] = self.scanDLL.sampleContrast(cSerialNumber)
            #p.ydataMean[p.ind] = self.scanDLL.sampleMean(cSerialNumber)
            #p.ydataSTD[p.ind] = self.scanDLL.sampleSTD(cSerialNumber)
          else:
            if self.state == "align":
              p.ydataContrast = p.ydataContrast[1:] + [self.scanDLL.dataContrast(cSerialNumber)]
              p.ydataMean = p.ydataMean[1:] + [self.scanDLL.dataMean(cSerialNumber)]
              p.ydataSTD = p.ydataSTD[1:] + [self.scanDLL.sampleSTD(cSerialNumber)]
        p.plotData(self.dataMode, self.nplot_updates_period)
    except: # pylint: disable=bare-except
      exc_type, exc_value, exc_traceback = sys.exc_info()# pylint: disable=redefined-outer-name
      traceback.print_tb(exc_traceback, limit=5, file=sys.stdout)
      self.log(repr(traceback.format_exception(exc_type, exc_value, exc_traceback)))

  def draw_plots(self):
    ''' Draw plots '''
    for p in self.plots:
      p.plotData(self.dataMode, 0)

  ########################################################
  #Prevent re-entering button click handlers on accidental double-clicks and queued clicks
  #self.buttonLock = False

  def set_buttonLock(self, lock):
    ''' re-entrance trigger '''
    self.buttonLock = lock
  #When parameter 'lock' is True, ";ock_button" locks ALL the buttons on the first entry
  #and returns True if called after buttons have been locked
  #To unlock, call this with False at the end of the click handler OR at the end of the thread created by the click handler

  def lock_buttons(self, lock):
    ''' Re-entrance handler '''
    if self.buttonLock and lock:
      print('Button locked')# An indicator of prevnted re-entrance. To be removed after completion of testing
      return True
    if lock:
      self.buttonLock = lock
    for b in self.buttons:
      b.blockSignals(lock)
      b.setEnabled(not lock)
      if not lock:
        QtCore.QTimer.singleShot(200, lambda: self.set_buttonLock(False))
    return False
  ########################################################

  def __del__(self):
    ''' ... '''
    if self.scanDLL is not None:
      self.scanDLL.close()

  def close(self):
    ''' Shut down the timer and cameras. '''
    if self.scanDLL is not None:
      try:
        #!!! self.scanDLL.stop()
        self.scanDLL.close()
        self.scanDLL = None
      except: # pylint: disable=bare-except
        exc_type, exc_value, exc_traceback = sys.exc_info()# pylint: disable=redefined-outer-name
        traceback.print_tb(exc_traceback, limit=5, file=sys.stdout)
        self.log(repr(traceback.format_exception(exc_type, exc_value, exc_traceback)))
    try:
      if self.thorLaserAvailable:
        self.thorLaserStatus = ctlThorITC.turnOffLaser(self.instrSeed,self.instrTA)
      if self.thorLaserAvailable:
        self.thorLaserStatus = ctlThorITC.turnOffTec(self.instrSeed,self.instrTA)
    except: # pylint: disable=bare-except
      if self.thorLaserAvailable:
        self.thorLaserStatus = ctlThorITC.closeLaserConn(self.instrSeed,self.instrTA,self.resMngr)
        if self.thorLaserStatus:
          self.thorLaserStatus = False # if closeLaserConn worked, set status to False
      self.log("Laser connection ")
    super().close()

  def closeEvent(self, event):
    ''' ... '''
    self.close()
    event.accept()

  def log(self, s):
    ''' Send a log message to the log text widget. '''
    self.logText.appendPlainText(s)
    self.logText.ensureCursorVisible()
    QtCore.QCoreApplication.instance().processEvents()  # this seems to work
    print(s)  # Print as well, for errors that cause crashes.
    self.logFile.write(s + '\n')

  def setButtonText(self, button, text):
    ''' ... '''
    button.setText(text)
    button.update()
    QtCore.QCoreApplication.instance().processEvents()

  def onExitButtonClick(self):
    ''' Multi-function button handler '''
    try:
      if self.thorLaserAvailable and self.thorLaserStatus:
        self.thorLaserStatus = ctlThorITC.turnOffLaser(self.instrSeed,self.instrTA)
      if self.state == "align":
        self.state = "stopping"
        self.exitButton.setEnabled(False)
        self.setButtonText(self.exitButton, "Stopping")
      else:
        if self.state == "calibrate":
          self.state = "stopping"
          self.exitButton.setEnabled(False)
          self.setButtonText(self.exitButton, "Cancelling")
        else:
          if self.state == "scan":
            self.state = "stopping"
            self.exitButton.setEnabled(False)
            self.setButtonText(self.exitButton, "Stopping")
          else:
            self.close()
    except: # pylint: disable=bare-except
      exc_type, exc_value, exc_traceback = sys.exc_info()# pylint: disable=redefined-outer-name
      traceback.print_tb(exc_traceback, limit=5, file=sys.stdout)
      self.log(repr(traceback.format_exception(exc_type, exc_value, exc_traceback)))

  # Data mode click handlers
  def onModeContrastClick(self):
    ''' ... '''
    self.dataMode = "Contrast"
    self.button_contrast.setDown(True)
    self.button_mean.setDown(False)
    self.button_std.setDown(False)
    self.draw_plots()
    QtCore.QCoreApplication.instance().processEvents()

  def onModeMeanClick(self):
    ''' ... '''
    self.dataMode = "Mean"
    self.button_contrast.setDown(False)
    self.button_mean.setDown(True)
    self.button_std.setDown(False)
    self.draw_plots()
    QtCore.QCoreApplication.instance().processEvents()

  def onModeSTDClick(self):
    ''' ... '''
    self.dataMode = "STD"
    self.button_contrast.setDown(False)
    self.button_mean.setDown(False)
    self.button_std.setDown(True)
    self.draw_plots()
    QtCore.QCoreApplication.instance().processEvents()

  def setModeButtons(self):
    ''' ... '''
    self.button_contrast.setDown(self.dataMode == "Contrast")
    self.button_mean.setDown(self.dataMode == "Mean")
    self.button_std.setDown(self.dataMode == "STD")
    QtCore.QCoreApplication.instance().processEvents()

  def onAlignClick(self):
    ''' Run the camera aligner. '''
    try:
      if self.lock_buttons(True):
        return False
      if self.thorLaserAvailable and self.thorLaserStatus:
        self.thorLaserStatus = ctlThorITC.turnOnLaser(self.instrSeed,self.instrTA)
      if self.thorLaserAvailable and not self.thorLaserStatus:
        self.log('ERROR: Thorlabs laser could not be turned on, connection disabled.')
        print('\a')
        return False
      self.state = "align"
      self.exitButton.setEnabled(False)
      # Create and write scan metadata, so the aligner can init the plugin.
      self.filename_ = 'align'
      scanDict = self.createScanDict()
      scanDict['hardwareParameters']['hwTrigger'] = False  # use sw triggering for aligner
      scanDict['cameraParameters']['numImages'] = 2500  # Run this longer between retriggers.
      scanDataDir = scanDict['fileParameters']['localScanDataDir']
      metadataFileName = scanDataDir + '/align_metadata.json'
      self.writeScanDict(scanDict, metadataFileName)

      # Point the log to the data output directory for the duration of the scan.
      savedLogFile = self.logFile
      self.logFile = open(scanDataDir + '/log.txt', 'a')  # closed on assigning to self.logFile

      # Initialize the scanner with metadata (creates CameraManager, OctopusManager, etc).
      cFilename = c_char_p(bytes(metadataFileName, 'utf-8'))
      if not self.scanDLL.init(cFilename, self.logFunc):
        self.log("ERROR: Scanner init failed.")
        self.logFile = savedLogFile
        return False

      # Create and run the aligner.
      self.alignerObj = aligner.Aligner(self, scanDict)
      if not self.alignerObj.good:
        self.logFile = savedLogFile
        return False
      if self.alignWindowModeless:
        self.alignerObj.show()
      else:
        self.alignerObj.exec_()  # return value ignored
      self.logFile = savedLogFile
    except: # pylint: disable=bare-except
      exc_type, exc_value, exc_traceback = sys.exc_info()# pylint: disable=redefined-outer-name
      traceback.print_tb(exc_traceback, limit=5, file=sys.stdout)
      self.log(repr(traceback.format_exception(exc_type, exc_value, exc_traceback)))
      self.log('Align exception.')
    finally:
      self.lock_buttons(False)
      self.state = "idle"
      self.exitButton.setEnabled(True)
      self.setButtonText(self.exitButton, "Exit")
      self.log('Aligner exited.')
      if self.thorLaserAvailable and self.thorLaserStatus:
        self.thorLaserStatus = ctlThorITC.turnOffLaser(self.instrSeed,self.instrTA)
    return True

  def onCalibrateClick(self):
    ''' Run the camera calibration procedure. '''
    # pylint: disable=too-many-return-statements
    try:
      if self.lock_buttons(True):
        return False
      if self.thorLaserAvailable and self.thorLaserStatus:
        self.thorLaserStatus = ctlThorITC.turnOnLaser(self.instrSeed,self.instrTA)
      if self.thorLaserAvailable and not self.thorLaserStatus:
        self.log('ERROR: Thorlabs laser could not be turned on, connection disabled.')
        print('\a')
        return False
      self.state = "calibrate"
      self.setButtonText(self.exitButton, "Cancel")
      self.exitButton.update()
      QtCore.QCoreApplication.instance().processEvents()
      self.log('Calibrating cameras ...')
      qtApp = QtCore.QCoreApplication.instance()

      # Create and write metadata.
      self.filename_ = 'calibrate'
      scanDict = self.createScanDict()
      scanDataDir = scanDict['fileParameters']['localScanDataDir']
      metadataFileName = scanDataDir + '/calibrate_metadata.json'
      self.writeScanDict(scanDict, metadataFileName)

      # Point the log to the data output directory for the duration of the scan.
      savedLogFile = self.logFile
      self.logFile = open(scanDataDir + '/log.txt', 'a')  # closed on assigning to self.logFile

      # Initialize the scanner with metadata (creates CameraManager, OctopusManager, etc).
      cFilename = c_char_p(bytes(metadataFileName, 'utf-8'))
      if not self.scanDLL.init(cFilename, self.logFunc):
        self.log("ERROR: Scanner init failed.")
        self.logFile = savedLogFile
        return False

      # Start the cameras. This clears the voxel buffers.
      if not self.scanDLL.start(c_bool(False)):
        self.log("ERROR: Can't start cameras.")
        self.logFile = savedLogFile
        return False

      # Collect frames for calibration.
      if not self.triggerAndWait():  # Error msgs will come from here.
        self.logFile = savedLogFile
        return False
      self.log('Collecting frames ...')
      self.scanDLL.beginAccumulation()  # This accumulates for all cameras.
      nsleep = 8
      gotFrames = False
      for _ in range(nsleep):
        time.sleep(1)
        #cancel calibration
        if self.state == "stopping":
          self.scanDLL.stop()
          self.log("Calibration cancelled by user")
          return False
        if self.scanDLL.anyBadFrames():
          self.log('ERROR: Bad frames detected. Try again.')
          self.scanDLL.stop()
          self.scanDLL.clearAccumulation()
          self.logFile = savedLogFile
          return False
        gotFrames = self.scanDLL.framesAccumulated()
        if gotFrames:
          pass  # Don't break; stopping early screws up data capture.
        qtApp.processEvents()

      # Beep, prompt, and stop the cameras.
      self.log('Frame averaging complete. Let up the trigger.')
      self.scanDLL.stop()
      time.sleep(1)

      # No frames? No good.
      if not gotFrames:
        self.log('ERROR: Frames not collected.')
        self.scanDLL.clearAccumulation()
        self.logFile = savedLogFile
        return False

      # Write the CSV.
      csvFile = scanDict['fileParameters']['localScanDataDir'] + "/data.csv"
      if self.writeCSV(csvFile):
        self.log('Wrote %s' % csvFile)
      else:
        self.log('WARNING: Error writing CSV file.')

      # Write the averaged frames. Without Gaussian flat-field correction, this is no longer
      # necessary, but we do this anyway to save the data for tracking camera issues.
      self.log('Writing frames ...')
      cFilenameBase = c_char_p(bytes(scanDataDir + '/cam', 'utf-8'))
      if not self.scanDLL.writeAverageFrames(cFilenameBase):
        self.log("WARNING: Writing of averaged frames failed.")
      self.scanDLL.clearAccumulation()

      # Proper parity is essential for calibration. Exception: No octopus -> no laser -> no light.
      dataGood = self.areDataGood(fixParity=False)
      if not dataGood:
        if self.noOctopusOk_:
          self.log('WARNING: Bad calibration data; ignoring (noOctopusOk).')
          self.calibrated_ = True
        else:
          self.log('ERROR: Bad calibration data. Try again.')
          self.calibrated_ = False
        print('\a')
      else:
        self.calibrated_ = True
        if windowPlatform:
          winsound.PlaySound('complete.wav', winsound.SND_ASYNC)

      self.calibrationFolder = ''
      if self.calibrated_:
        self.calibrationFolder = scanDataDir
      self.log('Done.')

      self.logFile = savedLogFile
    except: # pylint: disable=bare-except
      exc_type, exc_value, exc_traceback = sys.exc_info()# pylint: disable=redefined-outer-name
      traceback.print_tb(exc_traceback, limit=5, file=sys.stdout)
      self.log(repr(traceback.format_exception(exc_type, exc_value, exc_traceback)))
      self.log('Calibrate error.')
    finally:
      self.lock_buttons(False)
      self.state = "idle"
      self.exitButton.setEnabled(True)
      self.setButtonText(self.exitButton, "Exit")
      if self.thorLaserAvailable and self.thorLaserStatus:
        self.thorLaserStatus = ctlThorITC.turnOffLaser(self.instrSeed,self.instrTA)
    return True

  def onBackupClick(self):
    ''' Run a backup. '''
    try:
      if self.lock_buttons(True):
        return
      self.log('Starting backup of data with non-empty SubjectIDs.')
      src = outputDirectory()  # includes socket.gethostname()
      dst = "s3://owi-scan-data/" + socket.gethostname()
      subprocess.run(["aws", "s3", "sync", src, dst, '--exclude', '*', '--include', '*_scan_*', \
      '--include', '*_calibrate_*', '--include', '*_align_*'], check=False) # only backs up fold with non-empty subjectIDs
      self.log('Backup complete.')

    except: # pylint: disable=bare-except
      exc_type, exc_value, exc_traceback = sys.exc_info()# pylint: disable=redefined-outer-name
      traceback.print_tb(exc_traceback, limit=5, file=sys.stdout)
      self.log(repr(traceback.format_exception(exc_type, exc_value, exc_traceback)))
      self.log('Backup error.')
    finally:
      self.lock_buttons(False)

  def writeScanDict(self, scanDict, metadataFileName):
    ''' Write the scan metadata to disk. '''
    scanDataDir = scanDict['fileParameters']['localScanDataDir']
    try:
      pathlib.Path(scanDataDir).mkdir(parents=True, exist_ok=True)
      with open(metadataFileName, 'w') as f_local:
        json.dump(scanDict, f_local)
    except OSError as err:
      self.log("OS ERROR: {0}".format(err))
    except: # pylint: disable=bare-except
      exc_type, exc_value, exc_traceback = sys.exc_info()  # pylint: disable=redefined-outer-name
      traceback.print_tb(exc_traceback, limit=5, file=sys.stdout)
      self.log(repr(traceback.format_exception(exc_type, exc_value, exc_traceback)))
      self.log("Unexpected ERROR: %s" % sys.exc_info()[0])
      raise

  def triggerAndWait(self):
    ''' Trigger for hw (button) or sw (testing) data collection. '''
    ret = True
    if self.hwTriggered_:
      self.log("System ready for trigger")
      if not self.scanDLL.waitDataCollection(c_int(self.timeout_s_)):
        self.log('ERROR: Wait for data collection failed.')
        ret = False
    else:
      if self.scanDLL.triggerDataCollection():
        self.log('Sotware-triggered data collection.')
        ret = True
      else:
        self.log('ERROR: Software-triggered data collection failed.')
        ret = False
    return ret

  def showDiffuserConfirmation(self):
    ''' ... '''
    messagebox = QMessageBox(QMessageBox.Warning, " ", "", buttons = QMessageBox.Cancel | QMessageBox.Ok, parent=self)
    messagebox.setDefaultButton(QMessageBox.Cancel)
    messagebox.setStyleSheet("QPushButton {color:red; font-family: Arial; font-size:36px;}")
    messagebox.setIconPixmap(QPixmap('rsrc/FiberDiffuserCheck.jpg'))
    ret = messagebox.exec_()
    return ret == QMessageBox.Ok

  def onScanClick(self, _):# pylint: disable=too-many-return-statements
    ''' Start a scan. '''
    try:
      if self.lock_buttons(True):
        return False
      self.setModeButtons()
      if not self.calibrated_:
        if self.skip_calibration_for_debug:
          self.log("!!!!!!!!!!!!!!!!!!!!! Debug mode: Calibration skipped")
        else:
          self.log('The cameras are not calibrated yet. Run Calibrate first.')
          return False


      dllClassName = self.scanDLL.__class__.__name__
      if dllClassName == "CDLL":
        if not self.showDiffuserConfirmation():
          self.lock_buttons(False)
          self.state = "idle"
          return False

      qtApp = QtCore.QCoreApplication.instance()

      if self.thorLaserAvailable and self.thorLaserStatus:
        self.thorLaserStatus = ctlThorITC.turnOnLaser(self.instrSeed,self.instrTA)
      if self.thorLaserAvailable and not self.thorLaserStatus:
        self.log('ERROR: Thorlabs laser could not be turned on, connection disabled.')
        print('\a')
        return False
      self.log(self.calibrationFolder)

      self.state = "scan"
      self.setButtonText(self.exitButton, "Stop")
      self.exitButton.update()
      QtCore.QCoreApplication.instance().processEvents()
      self.log('Calibrating cameras ...')

      # Create and write scan metadata.
      self.filename_ = 'scan'
      scanDict = self.createScanDict()
      scanDataDir = scanDict['fileParameters']['localScanDataDir']
      metadataFileName = scanDataDir + '/scan_metadata.json'
      self.writeScanDict(scanDict, metadataFileName)

      # Point the log to the data output directory for the duration of the scan.
      savedLogFile = self.logFile
      self.logFile = open(scanDataDir + '/log.txt', 'a')  # closed on assigning to self.logFile

      # Log the scan start time.
      self.log('Scan start: %s' % time.asctime())

      # Initialize the scanner with metadata (creates CameraManager, OctopusManager, etc).
      cFilename = c_char_p(bytes(metadataFileName, 'utf-8'))
      if not self.scanDLL.init(cFilename, self.logFunc):  # This is a re-init, given aligner has run.
        self.log("ERROR: Scanner init failed.")
        self.logFile = savedLogFile
        return False

      # Iterate through the scan locations, prompting with new graphics.
      bloodflows = []  # bloodflows, L back to front, R front to back
      scanSuccessful = True
      for i in range(1, len(self.locations_) - 1):
        self.log('Capture at point %d ...' % i)
        #self.picLeft.setPixmap(QtGui.QPixmap(self.locations_[i][0]))
        #self.picRight.setPixmap(QtGui.QPixmap(self.locations_[i][1]))

        img = QtGui.QPixmap(self.locations_[i][0])
        img_height = img.height()
        self.picLeft.setPixmap(img.scaledToHeight(img_height / self.pixmap_scale))

        img = QtGui.QPixmap(self.locations_[i][1])
        img_height = img.height()
        self.picRight.setPixmap(img.scaledToHeight(img_height / self.pixmap_scale))

        qtApp.processEvents()  # Force a screen update.

        # N.B.: Subtle issue here: holding down the button while we loop around can trigger the next
        # location. So leave some time for the user to let the button up.

        # Loop over retries, if there are bad frames detected.
        locationSuccessful = False
        retries = 4
        for attempt in range(retries):
          self.log('Attempt %d of %d ...' % (attempt + 1, retries))
          # Start the cameras. This clears the voxel buffers.
          if not self.scanDLL.start(c_bool(False)):
            self.log("ERROR: Can't start cameras.")
            break
          # Trigger for hw (button) or sw (testing) data collection.
          if not self.triggerAndWait():  # error msgs will come from here
            break
          # Save a pair of dark/light test frames for each camera.
          if self.saveTestFrames_:
            for nth in range(self.scanDLL.rcam_NumCameras()):
              serialNumber = self.scanDLL.rcam_serialNumber(c_int(nth))
              filenameBase = scanDataDir + ("/loc%d_cam%d_" % (i, serialNumber))
              cFilenameBase = c_char_p(bytes(filenameBase, 'utf-8'))
              self.scanDLL.rcam_saveFrame(c_int(serialNumber), cFilenameBase, c_int(2))
          # Let the data capture process run. Periodically update the display.
          nSleep = int(scanDict['cameraParameters']['numImages'] * 2 /
                       scanDict['cameraParameters']['frameAcquisitionRate_Hz']) + 1
          self.log('  (capturing for %d secs)' % nSleep)

          self.clearPlots()
          et = QtCore.QElapsedTimer()
          et.start()
          while True:
            # if et.elapsed() > 8000:
            if et.elapsed() > (1000 * 1.1 * plot_array_size / scanDict['cameraParameters']['frameAcquisitionRate_Hz']):
              break
            if self.scanDLL is None:
              break
            self.update_plot()
            if not self.testAcceptBadFrames:
              if self.scanDLL.anyBadFrames():
                break
          et = None
          # Stop the cameras.
          if self.scanDLL is None:
            break
          if not self.scanDLL.stop():
            self.log("ERROR: Can't stop cameras.")
            break
          if self.state == "stopping":
            break
          # Wait 3 s for user to let go the button.
          self.log('Capture complete. Let up the trigger.')
          time.sleep(1)
          # If we got any bad frames, all bets are off w/dark,light frame tracking. Retry.
          if self.scanDLL.anyBadFrames():
            self.log('WARNING: Bad frames detected. Scan the same location again.')
            print('\a')
            time.sleep(2)
            continue
          # Run the data quality check.
          if not self.areDataGood(fixParity=True):
            self.log('WARNING: Bad scan data. Scan the same location again.')
            print('\a')
            time.sleep(2)
            continue
          locationSuccessful = True
          break
        # Failed in all retries? Outta here.
        if not locationSuccessful:
          scanSuccessful = False
          self.log('WARNING: Location %d failed after %d retries; continuing scan.' % (i, retries))
        # Write the CSV.
        if windowPlatform:
          winsound.PlaySound('complete.wav', winsound.SND_ASYNC)
        time.sleep(2)
        csvFile = scanDataDir + ("/location_%d.csv" % i)
        if not self.writeCSV(csvFile):
          break
        # Fit the data and print the result. (TODO(jfs): Display diff of right & left sides.)
        if locationSuccessful:
          alphaDb = self.fitBloodflow()
          self.log("Bloodflow: %g" % alphaDb)
          bloodflows.append(alphaDb)
        else:
          bloodflows.append(0.)
        #check if user clicked "Stop" button
        if self.state == "stopping":
          #self.scanDLL.stop()
          self.log("Scan stopped by user")
          return False
      # Stop the scanner (if it wasn't stopped above, after each location).
      self.scanDLL.stop()
      if scanSuccessful or i >= (len(self.locations_) - 2):
        img = QtGui.QPixmap(self.locations_[len(self.locations_) - 1][1])
        img_height = img.height()
        self.picRight.setPixmap(img.scaledToHeight(img_height / self.pixmap_scale))
        # Diff the lateral bloodflows. (ToDo(jfs): Put labels onscreen?)
        diffs = []
        for i in range(int(len(bloodflows) / 2)):
          diffs.append(abs(bloodflows[i] - bloodflows[len(bloodflows) - 1 - i]))
        self.log('Bloodflow differences: %s' % str(diffs))
        self.log('Scan end: %s' % time.asctime())
        self.logFile = savedLogFile
    except: # pylint: disable=bare-except
      exc_type, exc_value, exc_traceback = sys.exc_info()# pylint: disable=redefined-outer-name
      traceback.print_tb(exc_traceback, limit=5, file=sys.stdout)
      self.log(repr(traceback.format_exception(exc_type, exc_value, exc_traceback)))
      self.log('Scan error.')
    finally:
      self.lock_buttons(False)
      self.exitButton.setEnabled(True)
      self.setButtonText(self.exitButton, "Exit")
      if self.thorLaserAvailable and self.thorLaserStatus:
        self.thorLaserStatus = ctlThorITC.turnOffLaser(self.instrSeed,self.instrTA)
    return scanSuccessful

  def writeCSV(self, csvFile):
    ''' Write the captured data to a CSV file. '''
    if not self.scanDLL.writeCSV(bytes(csvFile, 'utf-8')):
      self.log("ERROR: Can't write %s." % csvFile)
      return False
    return True

  def areDataGood(self, fixParity):
    ''' Test for good scan data. '''
    for nth in range(self.scanDLL.rcam_NumCameras()):
      serialNumber = self.scanDLL.rcam_serialNumber(c_int(nth))
      if not self.scanDLL.checkData(c_int(serialNumber), c_bool(fixParity)):
        self.log('ERROR: Data check for camera %d failed.' % serialNumber)
        return False
    return True

  def fitBloodflow(self):
    ''' Return bloodflow fit to data currently residing in the plugin. '''
    if not self.scanDLL:
      return 0

    # Set up params for model fitting, one average per camera. Currently, this is single exposure.
    # ToDo(jfs): Simplify the interface to speckle_fcns for single exposure, if that's all we need.)
    mean, contrast, rhos = np.array([]), np.array([]), np.array([])
    for nth in range(self.scanDLL.rcam_NumCameras()):
      cameraID = self.scanDLL.rcam_serialNumber(c_int(nth))
      if str(cameraID) not in self._cameras:
        self.log('Camera %d not in camera info json' % cameraID)
        return 0
      info = self._cameras[str(cameraID)]
      rho = info['separation_mm']  # src-detector separation (mm)
      m = self.scanDLL.dataMean(c_int(cameraID))
      if m == 0:  # skip missing camera
        self.log('WARNING: Missing camera data (%d)? Data mean is 0.' % cameraID)
        return 0
      c = self.scanDLL.dataContrast(c_int(cameraID))
      mean = np.append(mean, m)
      contrast = np.append(contrast, c)
      rhos = np.append(rhos, rho)
      self.log('cam %d: mean %g, contrast %g' % (cameraID, m, c))
    if len(mean) == 0:  # test case: no cameras
      return 0

    # parameters
    exps = np.array([plot_array_size])*1e-6  # pulse widths / exposures (ToDo(jfs): Update this Gabor default.)
    params = {
      'n': 1.4,      # index of refraction
      'wv': 8.5e-4,  # wavelength, mm
      'mua': 0.005,  # absorption coefficient, mm^-1
      'musp': 1,     # reduced scattering coefficient, mm^-1
      'Db': 1e-6,    # alpha * diffusion coefficient of the scatterers, mm^2/s
      'dV2': 1,      # mean squared velocity of the scatterers, (mm/s)^2
      'Tmax': 1e-3,  # maximum exposure time being simulated, s
      'dtau': 1e-7,  # discretization of time step length, s
      'rho': rhos,   # distance between source and detector fibers, mm
      'exp': exps,   # exposure times, s
      'zb': 0.1,     # extrapolated boundary (for method of images), mm
      'beta': 0.25,  # constant from 0 to 1 depending on pixel size, polarization, and abberations, unitless
      'initial_guess': [40, 0.1], # for fitting I0 and mu_eff
      'usemueff': True,     # use the fitted mueff value in the flow fitting
      'fitbeta': 0,         # 1=fit for the beta parameter, 2=fit for beta and an offset
    }

    # arrange data in form expected by fitting functions (time, distance, pulse width)
    mean = np.expand_dims(np.expand_dims(mean, axis=1), axis=0)
    contrast = np.expand_dims(np.expand_dims(contrast, axis=1), axis=0)
    p = np.squeeze(params['rho'])
    nt = 1  # number of time points
    musp = params['musp']
    initial_guess = params['initial_guess']

    _, mua = speckle_fcns.fitattenuation(mean, p, nt, len(exps), musp, initial_guess)  # _ = mu_eff
    alphaDb, _, _ = speckle_fcns.fitflow_multi(params, contrast, nt, mua)  # _, _ = beta, offset


    return alphaDb[0][0]

  #
  # Scan Metadata
  #

  @staticmethod
  def getGitInfo():
    '''Record information about software version, branch, and latest updates.'''
    gitInfo = {}
    branch = subprocess.check_output(
        ['git', 'rev-parse', '--abbrev-ref', 'HEAD'], encoding='UTF-8').strip()
    gitInfo['branch'] = branch
    gitInfo['lastUpdate'] = subprocess.check_output(
        ['git', 'log', '-n', '1', '--pretty=format:%ar', "--", branch], encoding = 'UTF-8').strip()
    gitInfo['lastCommit'] = subprocess.check_output(
        ['git', 'log', '-n', '1', '--pretty=format:%h %s', "--", branch], encoding = 'UTF-8').strip()
    gitInfo['uncommittedMods'] = subprocess.check_output(
        ['git', 'status', '-s', '-uno'], encoding = 'UTF-8').strip()
    try:
      tagVersion = subprocess.check_output(
          ['git', 'describe', '--tags'], encoding = 'UTF-8').strip()
      gitInfo['tagVersion'] = tagVersion
    except subprocess.CalledProcessError:
      gitInfo['tagVersion'] = 'none'
    return gitInfo

  def createScanDict(self):
    ''' Generate scan dict from widget values. '''

    # AOM parameters
    aomParameters = {}
    aomParameters['AOM4Volt_V'] = 0.260 # voltage set by Sam, 10-01-2020
    aomParameters['AOM4Freq_MHz'] = 100

    # Camera parameters
    cameraParameters = {}
    rowTime_us = 15.8    # "Conversion Time" in CJ code
    resolutionX0, resolutionY0 = 0, 0 # For full frame: 0, 0 (half: 0,520)
    resolutionX1, resolutionY1 = 2712, 2080 # For full frame: 2712, 2080 (half: 2712, 1560)
    cameraParameters['pixelSize_um'] = 2.0
    cameraParameters['rowTime_us'] = rowTime_us
    cameraParameters['resolutionY0'] = resolutionY0  # [start location, pixels] [rows]
    cameraParameters['resolutionX0'] = resolutionX0  # [start location, pixels] [columns]
    cameraParameters['resolutionY1'] = resolutionY1  # [end location, pixels] [rows]
    cameraParameters['resolutionX1'] = resolutionX1  # [end location, pixels] [columns]
    cameraParameters['numImages'] = self.numImages_  # Octopus gates for 2x this, for dark/light pairs
    cameraParameters['frameAcquisitionRate_Hz'] = camFrameRate  # current reliable max w/gumstick

    # Camera Parameters for fake chopped system
    overlapTime_ms = 1.2  # Fix overlap time to ensure all pulse widths captured
    startOfLastRow_us = rowTime_us * (resolutionY1 - resolutionY0)
    cameraParameters['startOfLastRow_us'] = startOfLastRow_us
    cameraParameters['overlapTime_ms'] = overlapTime_ms
    cameraParameters['exposureTime_ms'] = (startOfLastRow_us * 10**-3) + overlapTime_ms
    cameraParameters["numCameras"] = self.scanDLL.rcam_NumCameras()
    # cameraParameters['cameraInfo'] = self._cameras # modified so that only the scanner's info gets recorded
    cameraParameters['cameraInfo'] = {}
    with open('cameras.json', 'r') as json_file:
      json_dict = json.load(json_file)
    scannerCameras = list(json_dict)
    # for i in range(len(scannerCameras)): # BH removed for pylint test
    for i, _ in enumerate(scannerCameras):
      scanName = socket.gethostname()
      if scannerCameras[i] == scanName.lower():
        for j in range(0,5):
          # testStr = scannerCameras[i+j] + ': ' + str(json_dict[scannerCameras[i+j]])
          cameraParameters['cameraInfo'][scannerCameras[i+j]] = json_dict[scannerCameras[i+j]]
    cameraParameters['calibrationFolder'] = self.calibrationFolder

    # Debug Parameters
    debugParameters = { 'noOctopusOk': self.noOctopusOk_ }

    # Delay parameters
    delayParameters = {}
    triggerOffset_s = 0.01  # 10ms offset to ensure no triggering on clock
    delayParameters["pulseWidth_s"] = 0.000100  # For Moglabs, at least. ThorLabs will be different.
    delayParameters["triggerOffset_s"] = triggerOffset_s
    delayParameters['frameValidDelay_s'] = triggerOffset_s - 0.001  # 1ms before FSIN begins
    delayParameters['frameValidWidth_s'] = 0.0025  # 2.5ms width; probably overkill
    delayParameters['chEDelay_s'] = triggerOffset_s - overlapTime_ms * 10**-3 + 0.1 * 10**-3  # always start pulse 0.1ms after rows start to overlap

    # Hardware parameters
    hardwareParameters = {}
    hardwareParameters["hwTrigger"] = self.hwTriggered_

    # Laser parameters
    laserParameters = {}
    laserParameters['laser'] = 'Moglabs850_Proto1_Zigzag'
    laserParameters['pseudoPulsed'] = 1  # Is this a pseudo-pulsed system

    # File Parameters
    gitInfo = self.getGitInfo()
    subjectID = self.subjectID.textEdit.text()
    fileParameters = { 'filename': self.filename_, 'subjectID': subjectID }
    fileParameters['swBranch'] = gitInfo['branch']
    fileParameters['swLastUpdate'] = gitInfo['lastUpdate']
    fileParameters['swLastCommit'] = gitInfo['lastCommit']
    fileParameters['swUncommittedMods'] = gitInfo['uncommittedMods']
    fileParameters['swTagVersion'] = gitInfo['tagVersion']
    strftime = datetime.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    rootFolder = strftime + '_' + self.filename_
    if subjectID != '':
      rootFolder += '_' + subjectID
    fileParameters['rootFolder'] = rootFolder
    # fileParameters['localScanDataDir'] = '.'
    fileParameters['localScanDataDir'] = outputDirectory() + '/' + rootFolder

    # Sample parameters
    sampleParameters = {
      'experimentNotes': self.experimentNotes.value(),
    }

    return {
      'AOMParameters': aomParameters,
      'cameraParameters': cameraParameters,
      'debugParameters': debugParameters,
      'delayParameters': delayParameters,
      'hardwareParameters': hardwareParameters,
      'laserParameters': laserParameters,
      'fileParameters': fileParameters,
      'sampleParameters': sampleParameters,
    }

def outputDirectory():
  ''' Return the directory where this machine's output files go. '''
  return 'c:/data_scans/' + socket.gethostname()

#
# main
#

def CheckConfig():
  ''' Check the machine configuration. Return True if it's what we expect. '''
  # Windows 10 version
  key = winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE, r"SOFTWARE\Microsoft\Windows NT\CurrentVersion")
  vn = int(winreg.QueryValueEx(key, "ReleaseID")[0])
  if vn not in (1909, 2004, 2009):
    print('ERROR: Wrong windows version. Scanner needs to be validated for this version.')
    return False

  # Python version
  if platform.python_version() != '3.7.6':
    print('ERROR: Wrong python version. Scanner needs to be validated for this version.')
    return False

  return True




if __name__ == "__main__":
  if not CheckConfig(): sys.exit(1)

  # Handle high resolution displays:
  if hasattr(QtCore.Qt, 'AA_EnableHighDpiScaling'):
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
  if hasattr(QtCore.Qt, 'AA_UseHighDpiPixmaps'):
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)

  # Create the main app. This is a singleton in Qt, retrievable elsewhere w/QApplication.instance().
  qapp = QtWidgets.QApplication(sys.argv)

  # Our Scanner is a top-level ("main") window. The rest is boilerplate setup.
  scanner = Scanner(scanDLL)
  scanner.show()
  frame_geom = scanner.frameGeometry()
  center = qapp.desktop().availableGeometry().center()
  frame_geom.moveCenter(center)
  scanner.move(frame_geom.topLeft())
  scanner.activateWindow()
  scanner.raise_()

  # Run the event loop / UI.
  try:
    qapp.exec_()
  except: # pylint: disable=bare-except
    exc_type, exc_value, exc_traceback = sys.exc_info()
    traceback.print_tb(exc_traceback, limit=5, file=sys.stdout)
    scanner.log(repr(traceback.format_exception(exc_type, exc_value, exc_traceback)))
