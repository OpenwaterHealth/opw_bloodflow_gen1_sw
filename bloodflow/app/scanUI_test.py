''' unit test for scanner class '''

import os
import sys

from PyQt5 import QtCore
from matplotlib.backends.qt_compat import QtWidgets

import scanUI

# mockScanDLL
class MockScanDLL():
  def __init__(self):
    self.numCameras_ = 0
    self.meanMean_ = 100.0
    self.darkMean_ = 99.0
    self.lightMean_ = 107.0
    self.serialNumbers_ = [ 3 ]
  def init(self, filename, logFunc):
    return True
  def log(self, s):
    print(s)
  def start(self, streamMode):
    return True
  def stop(self):
    return True
  def close(self):
    pass
  def triggerDataCollection(self):
    return True
  def rcam_NumCameras(self):
    return self.numCameras_
  def fx3_NumDevices(self, pid):
    return 0
  def writeCSV(self, filename):
    return True
  def dataMean(self, serialNumber):
    return 12.2
  def dataContrast(self, serialNumber):
    return 2.4
  def meanMean(self, serialNumber):
    return self.meanMean_
  def darkMean(self, serialNumber):
    return self.darkMean_
  def lightMean(self, serialNumber):
    return self.lightMean_
  def rcam_serialNumber(self, nth):
    return self.serialNumbers_[nth.value]  # nth is a c_int
  def anyBadFrames(self):
    return False
  def checkData(self, serialNumber, fixParity):
    return False

def test_onScanClickFailsWhenNotCalibrated():
  qapp = QtWidgets.QApplication(sys.argv)
  S = scanUI.Scanner(MockScanDLL())
  assert S.onScanClick(0) == False

def test_onScanClickSucceedsWhenCalibrated():
  qapp = QtWidgets.QApplication(sys.argv)
  S = scanUI.Scanner(MockScanDLL())
  S.calibrated_ = True
  S.numImages_ = 10
  assert S.onScanClick(0) == True

def test_createScanDict():
  qapp = QtWidgets.QApplication(sys.argv)
  S = scanUI.Scanner(MockScanDLL())
  S.experimentNotes.textEdit.setPlainText('brew ha ha')
  S.subjectID.textEdit.setText('thx1138')
  d = S.createScanDict()  # calls getGitInfo and writeScanDict
  assert d['sampleParameters']['experimentNotes'] == 'brew ha ha'
  assert d['laserParameters']['pseudoPulsed'] == 1
  assert d['cameraParameters']['resolutionY1'] == 2080
  assert d['cameraParameters']['resolutionX1'] == 2712
  assert d['fileParameters']['subjectID'] == 'thx1138'
  assert d['fileParameters']['swBranch'] is not None

def test_cameraInfo():
  qapp = QtWidgets.QApplication(sys.argv)
  S = scanUI.Scanner(MockScanDLL())
  assert len(S._cameras) > 0
  for key, info in S._cameras.items():
    assert key != ''
    if 'separation_mm' in info:
      assert info['separation_mm'] > 0

def testAreDataGoodFailsWhenCheckDataFails():
  qapp = QtWidgets.QApplication(sys.argv)
  dll = MockScanDLL()
  S = scanUI.Scanner(dll)
  dll.numCameras_ = 1
  assert S.areDataGood(fixParity=False) == False
