//
// Wrap C++ bloodflow image capture and analysis code in simple C calls for Python's CDLL package
//

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "bloodflow/component/execnode.h"
#include "bloodflow/component/frame.h"
#include "bloodflow/component/fx3.h"
#include "bloodflow/component/octopus.h"
#include "bloodflow/component/rcam.h"
#include "bloodflow/component/stddev.h"
#include "bloodflow/component/time.h"

#include "CameraManager.h"
#include "OctopusManager.h"
#include "VoxelData.h"
#include "VoxelSave.h"

using json = nlohmann::json;

//
// Wrapper
//

CameraManager* CameraMgr = nullptr;
OctopusManager* OctopusMgr = nullptr;

// How many debugging messages to print? Level 1 is the default, only important stuff.
int DebugLevel = 1;

// If you cannot afford an OctopusManager, one will be appointed for you ...
// ToDo(jfs): Octopus talks to hw only via FX3. It should be possible to use a real Octopus, with
// a fake FX3 (as in octopus_test), to get testing of OctopusMgr and Octopus together.
class MirandaOctopusManager: public OctopusManager {
 public:
  MirandaOctopusManager(): OctopusManager(nullptr /*Octopus*/) {}
  void TriggerDataCollection() override {  // Fake triggering by putting cameras in stream mode.
    CameraMgr->startAllCameras(true);
  }
  bool EnableSystemChannels(bool channelON) override { return true; }
  bool EnableFrameValid(bool validHIGH) override { return true; }
  bool ConfigureSystemChannels() override { return true; }
};

// Initialize the scanner given json metadata. If we've already initialized, restart OctopusMgr,
// as this is a cue we're switching modes from Align/calibrate to Scan.
extern "C" __declspec(dllexport) bool init(const char* fileName, void (*log)(const char*)) {
  // Read the metadata file.
  std::ifstream scanMetadata(fileName);
  if (scanMetadata.fail()) {
    log((std::string("ERROR: Can't open ") + fileName).c_str());
    return false;
  }
  json systemParameters = json::parse(scanMetadata);
  scanMetadata.close();
  log((std::string("Opened ") + fileName).c_str());

  // Disable all channels on Octopus (BH)
  if (Octopus::NumOctopi() > 0) {
    if (DebugLevel >= 1) log("Found an octopus; resetting.");
    OctopusMgr = new OctopusManager(new Octopus(new FX3()));
    if (DebugLevel >= 1) log("Initialized new OctopusManager.");
    if (OctopusMgr->init(systemParameters)) {
    } else {
      log("ERROR: OctopusManager init failed.");
      OctopusMgr = nullptr;
      return false;
    }
    OctopusMgr->EnableFrameValid(false);
    if (DebugLevel >= 1) log("Frame valid disabled.");
    OctopusMgr->EnableSystemChannels(false);  // Note: this should stop cameras & gates and stuff
    if (DebugLevel >= 1) log("Channels disabled.");
    delete OctopusMgr;
    OctopusMgr = nullptr;
  }
  
  // Create and intialize a CameraManager.
  if (CameraMgr == nullptr) {
    CameraMgr = new CameraManager([](){ return new Rcam(new FX3()); });
    if (!CameraMgr) return false;
    CameraMgr->logWith(log);
    try {
      if (!CameraMgr->init(systemParameters)) {
        delete CameraMgr;
        CameraMgr = nullptr;
        return false;
      }
    } catch(nlohmann::detail::exception e) {
      log((std::string("ERROR in CameraManager::init: ") + e.what()).c_str());
      return false;
    }
    log("CameraManager initialized.");
  } else if (systemParameters["fileParameters"]["filename"].get<std::string>() == "calibrate") {
    log("Reinitializing ROIs for calibration.");
    CameraMgr->ReinitROIs();
  }

  // Enable this to use fake OctopusManager if there's not a real one.
  bool noOctopusOk = systemParameters.contains("debugParameters")
      && systemParameters["debugParameters"]["noOctopusOk"].get<bool>();

  if (Octopus::NumOctopi() > 0) {
    if (DebugLevel >= 1) log("Found an octopus; initializing.");
    OctopusMgr = new OctopusManager(new Octopus(new FX3()));
    try {
      if (OctopusMgr->init(systemParameters)) {
        OctopusMgr->EnableFrameValid(true);
      } else {
        log("ERROR: OctopusManager init failed.");
        OctopusMgr = nullptr;
        delete CameraMgr;
        CameraMgr = nullptr;
        return false;
      }
    } catch(nlohmann::detail::exception e) {
      log((std::string("ERROR in OctopusManager::init: ") + e.what()).c_str());
      return false;
    }
  } else if (noOctopusOk) {
    if (DebugLevel >= 1) log("No octopus; using fake.");
    OctopusMgr = new MirandaOctopusManager();
  } else {
    log("ERROR: No octopus.");
    return false;
  }

  return true;
}

extern "C" __declspec(dllexport) bool start(bool streamMode) {
  // in lieu of Octopus checking every operation (WriteReg, etc)
  if (OctopusMgr->GetOctopus() && Octopus::NumOctopi() == 0) {
    CameraMgr->log("ERROR: Octopus failure.");
    return false;
  }

  if (!CameraMgr->startAllCameras(streamMode)) return false;
  OctopusMgr->ConfigureSystemChannels();  // ToDo: Test return, but alas this never returns false.
  OctopusMgr->EnableSystemChannels(true);  // Note: this will start cameras

  return true;
}

extern "C" __declspec(dllexport) bool triggerDataCollection() {
  OctopusMgr->TriggerDataCollection();
  return true;
}

extern "C" __declspec(dllexport) bool anyBadFrames() {
  return CameraMgr->anyBadFrames();
}

extern "C" __declspec(dllexport) bool stop() {
  if (DebugLevel >= 2) printf("dll.stop()\n");
  if (OctopusMgr->GetOctopus() && Octopus::NumOctopi() == 0) {
    CameraMgr->log("ERROR: Octopus failure.");
    CameraMgr->stopAllCameras();
    return false;
  }
  if (DebugLevel >= 2) printf("disabling OctopusMgr\n");
  OctopusMgr->EnableSystemChannels(false);  // Note: this should stop cameras
  if (DebugLevel >= 2) printf("stopAllCameras\n");
  return CameraMgr->stopAllCameras();
}

extern "C" __declspec(dllexport) void showMask(int serialNumber, bool show) {
  CameraMgr->showMask(serialNumber, show);
}

extern "C" __declspec(dllexport) bool writeCSV(const char* fileName) {
  return CameraMgr->writeCSVFile(std::string(fileName));
}

extern "C" __declspec(dllexport) bool checkData(int cameraID, bool fixParity) {
  return CameraMgr->checkData(cameraID, fixParity);
}

extern "C" __declspec(dllexport) float dataMean(int cameraID) {
  return CameraMgr->dataMean(cameraID);
}

extern "C" __declspec(dllexport) float dataContrast(int cameraID) {
  return CameraMgr->dataContrast(cameraID);
}

extern "C" __declspec(dllexport) float dataSTD(int cameraID) {
  return CameraMgr->dataSTD(cameraID);
}

extern "C" __declspec(dllexport) double meanMean(int cameraID) {
  return CameraMgr->meanMean(cameraID);
}

extern "C" __declspec(dllexport) double darkMean(int cameraID) {
  return CameraMgr->darkMean(cameraID);
}

extern "C" __declspec(dllexport) double lightMean(int cameraID) {
  return CameraMgr->lightMean(cameraID);
}

extern "C" __declspec(dllexport) float sampleContrast(int cameraID) {
  return CameraMgr->sampleContrast(cameraID);
}

extern "C" __declspec(dllexport) float sampleMean(int cameraID) {
  return CameraMgr->sampleMean(cameraID);
}

extern "C" __declspec(dllexport) float sampleSTD(int cameraID) {
  return CameraMgr->sampleSTD(cameraID);
}

extern "C" __declspec(dllexport) float *arrContrast(int cameraID) {
  return CameraMgr->arrayContrast(cameraID);
}

extern "C" __declspec(dllexport) float* arrMean(int cameraID) {
  return CameraMgr->arrayMean(cameraID);
}

extern "C" __declspec(dllexport) float* arrSTD(int cameraID) {
  return CameraMgr->arraySTD(cameraID);
}

extern "C" __declspec(dllexport) void clearAllCameras() {
  CameraMgr->clearAllCameras();
}

extern "C" __declspec(dllexport) void beginAccumulation() {
  CameraMgr->beginAccumulation();
}

extern "C" __declspec(dllexport) void clearAccumulation() {
  CameraMgr->clearAccumulation();
}

extern "C" __declspec(dllexport) bool framesAccumulated() {
  return CameraMgr->framesAccumulated();
}

extern "C" __declspec(dllexport) bool waitDataCollection(int timeout_s) {
  return CameraMgr->waitDataCollection(timeout_s);
}

extern "C" __declspec(dllexport) bool writeAverageFrames(const char* fileNameBase) {
  for (int i = 0, n = Rcam::NumCameras(); i < n; i++) {
    int cameraID = CameraMgr->serialNumber(i);
    if (!CameraMgr->writeAverageFrames(cameraID, std::string(fileNameBase) + std::to_string(cameraID))) {
      return false;
    }
  }
  return true;
}

extern "C" __declspec(dllexport) bool close() {
  if (OctopusMgr) {
    OctopusMgr->EnableSystemChannels(false);  // Note: this should stop cameras
    OctopusMgr->EnableFrameValid(false);
  }
  delete CameraMgr;
  CameraMgr = nullptr;
  delete OctopusMgr;
  OctopusMgr = nullptr;
  return true;
}

//
// Component::ExecNode
//

extern "C" __declspec(dllexport) int execnode_test() {
  ExecNode e;
  return (int)e.IsExecDone();
}

//
// Component::frame
//

extern "C" __declspec(dllexport) int frame_test() {
  Frame f1, f2;
  f1.SetTimestamp();
  Component::SleepMs(10);
  f2.SetTimestamp();
  printf("t1: %lld; t2: %lld\n", f1.timestamp_ms_, f2.timestamp_ms_);
  return f2.timestamp_ms_ > f1.timestamp_ms_;
}

//
// Component::fx3
//

extern "C" __declspec(dllexport) int fx3_NumDevices(int pid) {
  return FX3::NumDevices(pid);
}

extern "C" __declspec(dllexport) int fx3_Open_test(int pid, int nth) {
  FX3 fx3;
  int result = fx3.Open(pid, nth);
  if (result == 0) {
    printf("SerialNumber: %d\n", fx3.SerialNumber());
  }
  return result;
}

//
// Component::octopus
//

extern "C" __declspec(dllexport) int octopus_test() {
  int result = 0;
  if (Octopus::NumOctopi() > 0) {
    Octopus octo(new FX3());
    result = octo.Open();
  }
  return result;
}

//
// plugin::OctopusManager
//

const char* const TestJSON = R"foo({
  "cameraParameters": {
    "numCameras": 1,
    "numImages": 1,
    "frameAcquisitionRate_Hz": 15.0
  },
  "delayParameters": {
    "triggerOffset_s": 0.01,
    "frameValidDelay_s": 0.009,
    "frameValidWidth_s": 0.0025
  },
  "hardwareParameters": {
    "hwTrigger": 0
  }
}
)foo";

extern "C" __declspec(dllexport) int OctopusManager_test() {
  std::stringstream scanMetadata(TestJSON);
  json systemParameters = json::parse(scanMetadata);
  int result = 0;
  printf("NumOctopi: %d\n", Octopus::NumOctopi());
  if (Octopus::NumOctopi() > 0) {
    printf("Found an octopus.\n");
    OctopusManager omgr(new Octopus(new FX3()));
    bool initResult = omgr.init(systemParameters);
    if (!initResult) {
      printf("OctopusManager init failed\n");
    }
    result = initResult ? 0 : -1;
  }
  return result;
}

//
// Component::rcam
// No longer uses Rcam:: directly, but relies on CameraMgr being set up.
//

extern "C" __declspec(dllexport) int rcam_NumCameras() {
  return Rcam::NumCameras();
}

extern "C" __declspec(dllexport) int rcam_serialNumber(int nth) {
  return CameraMgr->serialNumber(nth);
}

extern "C" __declspec(dllexport) const char* rcam_model(int serialNumber) {
  return CameraMgr->getCamera(serialNumber)->Model();
}

extern "C" __declspec(dllexport) int rcam_size(int serialNumber) {
  Rcam* cam = CameraMgr->getCamera(serialNumber);
  Frame* frame = cam->GetConfig();
  return frame ? (frame->width << 16) | frame->height : 0;
}

extern "C" __declspec(dllexport) int rcam_saveFrame(int serialNumber, char* buffer, int nFrames) {
  if (DebugLevel >= 2) printf("rcam_saveFrame(%d)\n", serialNumber);
  CameraMgr->saveImages(serialNumber, nFrames, buffer);
  if (DebugLevel >= 2) printf("  waiting ...\n");
  return CameraMgr->waitForFrames(nFrames) ? 0 : -1;
}

//
// Component::stddev
//

extern "C" __declspec(dllexport) int stddev_test() {
  if (Rcam::NumCameras() == 0) return 0;

  Rcam cam(new FX3());
  int result = cam.Open(0);
  if (result != 0) return result;
  cam.Reset();
  cam.SetBLC(false);
  cam.SetExposure(0.001);
  cam.SetStream(true);
  Component::SleepMs(1000);

  StdDev stddev;
  stddev.AddProducer(&cam);

  cam.resize(30);
  cam.Start();

  // Let the pipeline run for a bit.
  // Output relies on a hacked stddev in lieu of VoxelSave.
  Component::SleepMs(2000);

  cam.Stop();

  while (!stddev.IsExecDone()) {
    Component::SleepMs(100);  // Wait for processing to complete before deleting any exec nodes
  }

  cam.Close();

  return 0;
}

//
// Component::time
//

extern "C" __declspec(dllexport) time_t SteadyClockTimeMs() {
  return Component::SteadyClockTimeMs();
}

extern "C" __declspec(dllexport) int SleepMs(int ms) {
  Component::SleepMs(ms);
  return 0;
}

//
// plugin::VoxelSave
//

extern "C" __declspec(dllexport) int voxelsave_test() {
  if (Rcam::NumCameras() == 0) return 0;

  Rcam cam(new FX3());
  int result = cam.Open(0);
  if (result != 0) return result;
  cam.Reset();
  cam.SetBLC(true);
  cam.SetExposure(0.001);
  cam.SetStream(true);
  Component::SleepMs(1000);

  StdDev stddev;
  stddev.AddProducer(&cam);

  VoxelSave vsave(
    [](const VoxelData& v){ printf("mean = %g, contrast = %g\n", v.imageMean, v.speckleContrast); },
    [](Frame*){}, 0.0);
  vsave.AddProducer(&stddev);

  cam.resize(30);
  cam.Start();
  Component::SleepMs(100);

  // Let the pipeline run for a bit.
  Component::SleepMs(2000);

  cam.Stop();

  while (!vsave.IsExecDone()) {
    Component::SleepMs(100);  // Wait for processing to complete before deleting any exec nodes
  }

  cam.Close();

  return 0;
}
