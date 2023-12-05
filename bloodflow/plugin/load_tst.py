import json
import os
import unittest

from ctypes import CDLL, WINFUNCTYPE, c_bool, c_char_p, c_int

logFuncType = WINFUNCTYPE(c_int, c_char_p)
def logCallback(s):
  print(str(s, 'utf-8') + '\n')
  return 0
logFunc = logFuncType(logCallback)  # for extra credit: lambda s: log(s)

class LoadTest(unittest.TestCase):
  def setUp(self):
    self.scanDLLPath = 'bloodflow/plugin/scan.dll'
    # self.scanDLLPath = '../../bazel-bin/bloodflow/plugin/scan.dll'  # for running outside bazel
    self.scanDLL = CDLL(self.scanDLLPath)
    self.scanDLL.init.restype = c_bool
    self.scanDLL.start.restype = c_bool
    self.scanDLL.rcam_model.restype = c_char_p
    self.scanDLL.writeCSV.restype = c_bool
    self.scanDLL.close.restype = c_bool
    print('Loaded scan.dll from', self.scanDLLPath)

  #
  # CameraManager / Wrapper
  #
  def testCameraManager(self):
    print('testCameraManager')
    scanDict = {
      "cameraParameters": { "exposureTime_ms": 1.0, "gain": 1.0, "cameraInfo": {} },
      "debugParameters": { "noOctopusOk": False },
      "delayParameters": { "pulseWidth_s": 0.000100 },
      "fileParameters": { "filename": "none" }
    }
    metadataFileName = b'scan_metadata.json'
    with open(metadataFileName, 'w') as f_local:  # closes at end of scope
      json.dump(scanDict, f_local)
    self.assertFalse(self.scanDLL.init(c_char_p(metadataFileName), logFunc))  # fails w/o octopus
    scanDict['debugParameters']['noOctopusOk'] = True
    with open(metadataFileName, 'w') as f_local:  # closes at end of scope
      json.dump(scanDict, f_local)
    self.assertTrue(self.scanDLL.init(c_char_p(metadataFileName), logFunc))
    self.assertTrue(self.scanDLL.start(c_bool(True)))  # streamMode for testing w/o trigger
    print('sleeping for 2 secs')
    self.scanDLL.SleepMs(c_int(2000))
    csvFileName = b'data.csv'
    self.assertTrue(self.scanDLL.writeCSV(csvFileName))
    self.assertTrue(self.scanDLL.close())
    os.unlink(metadataFileName)

  #
  # Component::ExecNode
  #
  def testExecNode(self):
    print('testExecNode')
    self.assertTrue(self.scanDLL.execnode_test() == 1)

  #
  # Component::frame
  #
  def testFrame(self):
    print('testFrame')
    self.assertTrue(self.scanDLL.frame_test())

  #
  # Component:fx3
  #
  def testFX3NumDevices(self):
    numDevices = self.scanDLL.fx3_NumDevices(c_int(0x4F10))
    print('fx3_NumDevices:', numDevices)
    self.assertTrue(numDevices >= 0)  # hard to make this anything else, w/o cameras or mockCypress
  def testFXOpen(self):
    numDevices = self.scanDLL.fx3_NumDevices(c_int(0x4F10))
    result = self.scanDLL.fx3_Open_test(c_int(0x4F10), c_int(0))
    print('fx3_Open_test:', result)
    self.assertTrue(result == (0 if numDevices > 0 else -1))

  #
  # Component::octopus
  #
  def testOctopus(self):
    print('testOctopus')
    self.assertTrue(self.scanDLL.octopus_test() == 0)

  #
  # plugin::OctopusManager
  #
  def testOctopusManager(self):
    print('testOctopusManager')
    self.assertTrue(self.scanDLL.OctopusManager_test() == 0)

  #
  # Component::rcam
  #
  def testRcam(self):
    numCameras = self.scanDLL.rcam_NumCameras()
    print('rcam_NumCameras:', numCameras)
    self.assertTrue(numCameras >= 0)
    # ToDo(jfs): Move these tests to testCameraManager().
    # for nth in range(numCameras):
    #   serialNumber = self.scanDLL.rcam_open(c_int(nth))
    #   self.assertTrue(serialNumber > 0)
    #   frame_size = self.scanDLL.rcam_size(c_int(serialNumber))
    #   width = (frame_size >> 16) & 0xFFFF
    #   height = frame_size & 0xFFFF
    #   print("Camera", nth, "opened: serialNumber:", serialNumber, ";", width, "x", height)
    #   self.assertTrue(width > 0 and height > 0)
    #   model = self.scanDLL.rcam_model(c_int(serialNumber))
    #   self.assertTrue(model is not None)
    #   print("  model:", model)
    #   self.assertTrue(self.scanDLL.rcam_start(c_int(serialNumber)) == 0)
    #   print('started')
    #   self.scanDLL.rcam_close(c_int(serialNumber))
    #   print('closed')

  #
  # Component::StdDev
  #
  def testStdDev(self):
    print('testStdDev')
    self.assertTrue(self.scanDLL.stddev_test() == 0)

  #
  # Component:time
  #
  def testSteadyClockTimeMs(self):
    t = self.scanDLL.SteadyClockTimeMs()
    print('SteadyClockTimeMs =', t)
    self.assertTrue(t != 0)
  def testSleepMs(self):
    print('testSleepMs')
    self.assertTrue(self.scanDLL.SleepMs(1000) == 0)

  #
  # plugin::VoxelSave
  #
  def testVoxelSave(self):
    print('testVoxelSave')
    self.assertTrue(self.scanDLL.voxelsave_test() == 0)

if __name__ == "__main__":
  path = os.getenv('PATH')
  for p in path.split(';'):
    pp = p + '/tiff.dll'
    print(os.path.exists(pp), pp)
  os.environ['PATH'] += ';./bloodflow/plugin'
  unittest.main()
