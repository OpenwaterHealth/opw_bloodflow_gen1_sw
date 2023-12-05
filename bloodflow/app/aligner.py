''' Camera aligner for bloodflow UI '''
import sys
#import traceback

from ctypes import c_bool, c_char_p, c_int
from matplotlib.backends.qt_compat import QtCore, QtGui, QtWidgets

class Aligner(QtWidgets.QDialog):
  ''' Camera alignment interface '''
  def __init__(self, scanner, scanDict):
    ''' Initialize, given scanner for DLL, etc. '''
    super().__init__(scanner)
    self._scanner = scanner
    self._scanDLL = scanner.scanDLL
    self._scanDict = scanDict
    self._cameras_layout = QtWidgets.QVBoxLayout()
    self._layout = QtWidgets.QVBoxLayout(self)
    self._cameras = []
    self._images = []
    self._means = []
    self._meanMeans = []
    self._darkMeans = []
    self._lightMeans = []
    self._contrasts = []
    self._timer = QtCore.QTimer(self)
    self._timer.timeout.connect(self.updateCameraImage)
    self._zoom = 4
    self._imageBuffer = None
    self._periodMs = 500  # was 1000, w/tiff temp files
    self._useDataPlots = False
    #Scale widgets relative to the base screen size
    screen_sz = QtCore.QCoreApplication.instance().primaryScreen().size()
    base_screen_h = 1260
    pixmap_scale = base_screen_h / screen_sz.height()
    self._zoom *= pixmap_scale
    self.good = self.setupCameras()  # calls start()
    if not self.good:
      self._timer.stop()
      self._scanDLL.stop()
    self._layout.addLayout(self._cameras_layout)
    stop_button = self.StopButton()
    self._layout.addWidget(stop_button)
    if self._useDataPlots:
      self._scanner.createDataPlotForCameras()
      self._scanner.clearPlots()

  def centerWindow(self):
    ''' ... '''
    frame_geom = self.frameGeometry()
    qapp = QtWidgets.QApplication(sys.argv)
    center = qapp.desktop().availableGeometry().center()
    frame_geom.moveCenter(center)
    self.move(frame_geom.topLeft())

  def cleanClose(self):
    ''' ... '''
    self._timer.stop()
    #!!!
    self._timer = None
    if self._scanDLL is not None:
      self._scanDLL.stop()

  def closeEvent(self, event):
    ''' ... '''
    self.cleanClose()
    event.accept()

  def onStopClick(self):
    ''' ... '''
    super().close()

  def StopButton(self):
    ''' Multi-function Exit / Close / Cancel button '''
    button = QtWidgets.QPushButton("Stop")
    button.clicked.connect(self.onStopClick)
    button.setStyleSheet("background-color: lightgrey;")
    bf = button.font()
    bf.setPointSize(bf.pointSize()*4)
    button.setFont(bf)
    return button

  def log(self, s):
    ''' Convenience function: route log call to scanner '''
    self._scanner.log(s)

  def start(self):
    ''' Start the cameras. '''
    if not self._scanDLL.start(c_bool(False)):
      self.log("ERROR: Can't start cameras.")
      return False
    if not self._scanDLL.triggerDataCollection():
      self.log('ERROR: Data collection trigger failed.')
      return False
    return True

  def stop(self):
    ''' Stop the cameras. '''
    if not self._scanDLL.stop():
      self.log('ERROR: ScanDLL stop() failed.')
      return False
    self.saveCameraImages()
    return True

  def addField(self, title, vbox):
    ''' Add a field with a label. '''
    vbox.addWidget(QtWidgets.QLabel(title))  # add the title label
    label = QtWidgets.QLabel()  # add the value label
    vbox.addWidget(label)
    return label

  def getCameraPixmap(self, serialNumber):
    ''' Get a camera frame as a pixmap. '''
    cSerialNumber = c_int(serialNumber)
    cBuffer = c_char_p(self._imageBuffer.bits().__int__())
    saveResult = self._scanDLL.rcam_saveFrame(cSerialNumber, cBuffer, c_int(1))
    if saveResult != 0:
      self.log('ERROR %d saving frame for camera %d' % (saveResult, serialNumber))
      # Return an empty pixmap.
    height = self._imageBuffer.height()
    return QtGui.QPixmap.fromImage(self._imageBuffer.scaledToHeight(height / self._zoom))

  def cameraLayout(self, serialNumber):
    ''' Return a layout for one camera, adding to _cameras & other arrays. '''

    model = self._scanDLL.rcam_model(c_int(serialNumber))
    frame_size = self._scanDLL.rcam_size(c_int(serialNumber))
    width = (frame_size >> 16) & 0xFFFF
    height = frame_size & 0xFFFF
    if self._imageBuffer is None:  # ToDo(jfs): Check all are the same size.
      self._imageBuffer = QtGui.QImage(width, height, QtGui.QImage.Format_Grayscale8)
    self.log("Camera %d, model %s, size %d x %d." % (serialNumber, model, width, height))
    self._cameras.append(serialNumber)
    vbox = QtWidgets.QVBoxLayout()  # info/button
    vbox.addWidget(QtWidgets.QLabel(str(serialNumber)))  # Add a label for the serial number.
    vbox.addWidget(QtWidgets.QLabel('%d x %d' % (width, height)))

    self._means.append(self.addField('mean', vbox))
    self._contrasts.append(self.addField('contrast', vbox))
    self._meanMeans.append(self.addField('meanMean', vbox))
    self._darkMeans.append(self.addField('darkMean', vbox))
    self._lightMeans.append(self.addField('lightMean', vbox))

    checkbox = QtWidgets.QCheckBox()
    checkbox.setText('Show Mask')
    # In PyQt this callback takes the value as first arg, so add serialNumber as second.
    checkbox.clicked.connect(lambda val, sn=serialNumber: self.onShowMask(sn, val))
    vbox.addWidget(checkbox)

    label = QtWidgets.QLabel()  # Add a label for the image.
    self._images.append(label)
    pixmap = self.getCameraPixmap(serialNumber)
    if pixmap is None:
      return None
    label.setPixmap(pixmap)

    hbox = QtWidgets.QHBoxLayout()
    hbox.addLayout(vbox)
    hbox.addWidget(label)
    return hbox

  def onShowMask(self, serialNumber, val):
    ''' Turn the mask on or off. '''
    self._scanDLL.showMask(c_int(serialNumber), c_bool(val))

  def setupCameras(self):
    ''' Find the attached cameras and set them up. '''

    # Clear any existing image widgets.
    while self._cameras_layout.count():
      w = self._cameras_layout.itemAt(0).widget()
      w.setText('')
      self._cameras_layout.removeWidget(w)

    # Start the cameras here so we can grab images.
    if not self.start():
      return False

    # Grab an image from each camera. Create a list of cameras and image widgets during the process.
    self._cameras = []
    self._images = []
    self._means = []
    self._contrasts = []
    self._meanMeans = []
    self._darkMeans = []
    self._lightMeans = []
    numCameras = self._scanDLL.rcam_NumCameras()
    self._scanner.log('Found %d camera(s)' % numCameras)
    pair = QtWidgets.QHBoxLayout()  # two cameras per row
    for nth in range(numCameras):
      serialNumber = self._scanDLL.rcam_serialNumber(c_int(nth))
      if serialNumber <= 0:
        self.log("ERROR: Can't open camera, index %d" % nth)
        continue
      layout = self.cameraLayout(serialNumber)
      if layout is None:
        return False
      pair.addLayout(layout)
      if nth % 2 == 1 or nth == numCameras - 1:  # on odd or last camera ...
        self._cameras_layout.addLayout(pair)
        if nth < numCameras - 1:
          pair = QtWidgets.QHBoxLayout()
    self.adjustSize()

    # Kick off the update timer.
    self._timer.start(self._periodMs)

    return True

  def updateCameraImage(self):
    ''' Update the image for a camera. '''
    if self._useDataPlots:
      self._scanner.update_plot()
    # Grab a new image for each camera label widget.
    n = 0
    for serialNumber in self._cameras:
      cSerialNumber = c_int(serialNumber)
      dataMean = self._scanDLL.dataMean(cSerialNumber)
      if dataMean == -1:  # triggerDataCollection / numImages expired; restart
        self._timer.stop()
        self.log('retriggering')
        if (not self.stop()) or (not self.start()):  # Emits error message.
          self.log('ERROR: Aligner restart failed.')
          return  # W/o _timer.start, aligner is effectively off.
        self._timer.start(self._periodMs)
        return  # Get out cleanly, and wait for the next update.
      dataContrast = self._scanDLL.dataContrast(cSerialNumber)
      pixmap = self.getCameraPixmap(serialNumber)
      if pixmap is not None:
        self._images[n].setPixmap(pixmap)
        # Get the stats and update the labels.
        numder_of_digits = 3
        self._means[n].setText(self._means[n].locale().toString(dataMean, 'f', numder_of_digits))
        self._contrasts[n].setText(self._contrasts[n].locale().toString(dataContrast, 'f', numder_of_digits))
        self._meanMeans[n].setText(self._meanMeans[n].locale().toString(self._scanDLL.meanMean(cSerialNumber), 'f', numder_of_digits))
        self._darkMeans[n].setText(self._darkMeans[n].locale().toString(self._scanDLL.darkMean(cSerialNumber), 'f', numder_of_digits))
        self._lightMeans[n].setText(self._lightMeans[n].locale().toString(self._scanDLL.lightMean(cSerialNumber), 'f', numder_of_digits))
      n += 1
    self._scanDLL.clearAllCameras()  # Clear voxel accumulators.

  def saveCameraImages(self):
    ''' Save the latest image for each camera. '''
    scanDataDir = self._scanDict['fileParameters']['localScanDataDir']
    n = 0
    for serialNumber in self._cameras:
      filename = scanDataDir + ('/align_camera_%d.tiff' % serialNumber)
      if self._images[n].pixmap().save(filename):
        self.log('Saved %s' % filename)
      else:
        self.log('ERROR saving image for camera %d' % serialNumber)
      n += 1
