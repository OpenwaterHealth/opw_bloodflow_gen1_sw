''' unit test for aligner class '''

import os
import sys

from PyQt5 import QtCore
from matplotlib.backends.qt_compat import QtWidgets

import aligner

# mockScanner
class MockScanner(QtWidgets.QMainWindow):
  def __init__(self, scanDLL_):
    super().__init__()
    self.scanDLL = scanDLL_
    self.logFunc = None  # uses WINFUNCTYPE, etc in Scanner
  def log(self, s):
    print(s)

# mockScanDLL
class MockScanDLL():
  def __init__(self):
    self.numCameras_ = 0
    self.starts_ = 0
    self.saveFrames_ = 0
  def init(self, metadataFileName, log):
    return True
  def rcam_NumCameras(self):
    return self.numCameras_
  def rcam_serialNumber(self, nth):
    return 999
  def rcam_open(self, nth):
    return 301  # fake serialNumber
  def rcam_model(self, nth):
    return 'Hi, Max!'
  def rcam_size(self, nth):
    return 0x01000100
  def rcam_start(self, serialNumber):
    return True
  def rcam_saveFrame(self, serialNumber, cBuffer, nFrames):
    self.saveFrames_ += 1
    return 0
  def triggerDataCollection(self):
    return True
  def clearAllCameras(self):
    return True
  def dataMean(self, serialNumber):
    return 500.0
  def dataContrast(self, serialNumber):
    return 0.62
  def meanMean(self, serialNumber):
    return 3.0
  def darkMean(self, serialNumber):
    return 2.2
  def lightMean(self, serialNumber):
    return 3.8
  def start(self, streamMode):
    self.starts_ += 1
    return True
  def stop(self):
    pass
  def close(self):
    pass

def testCreateAndDestroyWithNoCameras():
  qapp = QtWidgets.QApplication(sys.argv)
  scanDLL = MockScanDLL()
  S = MockScanner(scanDLL)
  A = aligner.Aligner(S, {})
  assert True

def testCreateAndDestroyWithOneCamera():
  qapp = QtWidgets.QApplication(sys.argv)
  scanDLL = MockScanDLL()
  scanDLL.numCameras_ = 1
  S = MockScanner(scanDLL)
  A = aligner.Aligner(S, {})
  assert len(A._images) == 1

def testUpdate():
  qapp = QtWidgets.QApplication(sys.argv)
  scanDLL = MockScanDLL()
  scanDLL.numCameras_ = 1
  S = MockScanner(scanDLL)  # calls saveFrame the first time
  A = aligner.Aligner(S, {})
  A.updateCameraImage()  # should result in another call to saveFrame
  assert scanDLL.saveFrames_ == 2
