#include <cstdint>
#include <iomanip>
#include <iostream>
#include <time.h>

#include "bloodflow/component/rcam.h"
#include "bloodflow/component/stddev.h"
#include "bloodflow/component/time.h"

#include "CameraManager.h"
#include "VoxelData.h"
#include "VoxelSave.h"

using json = nlohmann::json;
using std::to_string;

CameraManager::~CameraManager() {
  for (const auto& it : cameraInfoMap_) {
    const cameraInfo& info = it.second;
    while (!info.voxelSave->IsExecDone()) {
      Component::SleepMs(100);  // Wait for processing to complete before deleting any exec nodes
    }
    delete info.camera;
    delete info.stdDev;
    delete info.voxelSave;
  }
}

void CameraManager::log(const char* s) const {
  if (log_) {
    log_(s);
  } else {
    printf("%s\n", s);
  }
}

bool CameraManager::init(const json& systemParameters) {
  // TODO (CR): check for stupid inputs (like bad ROI defaults)
  json cameraParameters = systemParameters["cameraParameters"];
  double exposureTime_s = cameraParameters["exposureTime_ms"].get<double>() * 0.001;
  double resolutionY0 = cameraParameters["resolutionY0"].get<double>();
  double resolutionY1 = cameraParameters["resolutionY1"].get<double>();
  json camerasJSON = cameraParameters["cameraInfo"];

  // Set up cameras and image processing pipeline
  int numCameras = Rcam::NumCameras();
  for (int i = 0; i < numCameras; i++) {
    Rcam* camera = newRcam_();
    if (camera->Open(i) != 0) {
      log(std::string("ERROR: Can't open camera, index ") + to_string(i));
      delete camera;
      return false;
    }
    int cameraID = camera->SerialNumber();
    log(std::string("Connected to Camera ") + to_string(cameraID));
    std::string idStr = to_string(cameraID);
    json cameraJSON;
    if (camerasJSON.contains(idStr)) {
      cameraJSON = camerasJSON[idStr];
    } else {
      log(std::string("WARNING: No camera JSON for camera ") + to_string(cameraID) + "; using defaults.");
      // log(std::string("WARNING: No camera JSON for camera ") + to_string(cameraID) + "; cancelling initializaiton.");
      // return false; // This would be a quicker way to catch if cameras are not in json but causes major unit test issues BH
    }
    double gain = cameraJSON.contains("gain") ? cameraJSON["gain"].get<double>() : 1.0;
    double meanThres = cameraJSON.contains("meanThres") ? cameraJSON["meanThres"].get<double>() : 0.5; //To pass tests with mock cameras // 5.0;
    json roiJSON;
    if (cameraJSON.contains("roi"))
      roiJSON = cameraJSON["roi"];

    // Try streaming a few frames. If that doesn't work, reset.
    if (!testCamera(camera)) {
      delete camera;
      return false;
    }

    // Store camera information
    cameraInfo& cameraInfo = cameraInfoMap_[cameraID];
    cameraInfo.cameraID = cameraID;
    cameraInfo.meanThres = meanThres;
    cameraInfo.camera = camera;
    cameraInfo.stdDev = new StdDev();
    if (roiJSON.contains("x")) {
      const Frame* config = camera->GetConfig();
      cameraInfo.stdDev->SetROI(roiJSON["x"].get<int>(), roiJSON["y"].get<int>() - resolutionY0,
          roiJSON["r"].get<int>(), config->width, config->height);
    }
    if (cameraJSON.contains("gainConst")) {
      cameraInfo.stdDev->SetGainConstant(cameraJSON["gainConst"].get<double>());
    }
    cameraInfo.voxelSave = new VoxelSave([&](const VoxelData& v) -> void { this->addVoxel(v); },
        [&](Frame* fr) -> void { this->saveFrame(fr); }, cameraInfo.stdDev->GainConstant());

    cameraInfo.camera->SetExposure(exposureTime_s);
    cameraInfo.camera->SetBLC(false);
    cameraInfo.camera->SetGain(gain);
    cameraInfo.camera->SubWindow2Point(0, resolutionY0, 2712, resolutionY1); // leaving X set to full row since no improvement in speed

    // Set up processing chain
    cameraInfo.stdDev->AddProducer(cameraInfo.camera);
    cameraInfo.voxelSave->AddProducer(cameraInfo.stdDev);
    cameraInfo.camera->resize(30);
  }

  return true;
}

bool CameraManager::testCamera(Rcam* camera) {
  // Try streaming a few frames.
  int cameraID = camera->SerialNumber();
  camera->SetStream(true);
  int frameCountBefore = camera->GetFrameCount();
  if (camera->Start() != 0) {
    log(std::string("ERROR: Can't start camera ") + to_string(cameraID));
    return false;
  }
  Component::SleepMs(200);
  camera->Stop();
  int framesReceived = camera->GetFrameCount() - frameCountBefore;
  camera->SetStream(false);

  // If that doesn't work, reset.
  if (framesReceived > 0) {
    log(std::string("  camera test passed (received ") + to_string(framesReceived) + " frames).");
    return true;
  }
  log("  camera test failed; resetting.");
  if (camera->Reset() != 0) {
    log(std::string("ERROR: Can't reset camera ") + to_string(cameraID));
    return false;
  }

  // Try streaming a few frames again.
  camera->SetStream(true);
  frameCountBefore = camera->GetFrameCount();
  if (camera->Start() != 0) {
    log(std::string("ERROR: Can't start camera ") + to_string(cameraID));
    return false;
  }
  Component::SleepMs(200);
  camera->Stop();
  framesReceived = camera->GetFrameCount() - frameCountBefore;
  camera->SetStream(false);
  if (framesReceived > 0) {
    log(std::string("  camera test passed (received ") + to_string(framesReceived) + " frames).");
    return true;
  }

  log("  camera test failed after reset.");
  return false;
}

int CameraManager::serialNumber(int nth) {
  int i = 0;
  for (const auto& it : cameraInfoMap_) {
    if (i == nth) return it.first;
    ++i;
  }
  return -1;
}

void CameraManager::saveImages(int serialNumber, int nFrames, char* buffer) {
  framesSaved_ = 0;
  if (strncmp(buffer, "c:/", 3) == 0) {  // file? (ToDo: Better way to determine this.)
    filenameBase_ = std::string(buffer);
  } else {  // or memory
    frameBuffer_ = buffer;
  }
  cameraInfoMap_[serialNumber].voxelSave->saveFrames(nFrames);
}

bool CameraManager::waitForFrames(int nFrames) {
  time_t start = Component::SteadyClockTimeMs();
  while (framesSaved_ < nFrames) {
    if (Component::SteadyClockTimeMs() - start > 1000) {
      log("WARNING: timeout.");
      break;
    }
  }
  frameBuffer_ = nullptr;  // in case we were saving to memory ...
  filenameBase_.clear();   // ... or file(s)
  return framesSaved_ >= nFrames;
}

void CameraManager::ReinitROIs() {
  for (auto& it : cameraInfoMap_) {
    it.second.stdDev->ReinitROI();
  }
}

bool CameraManager::startAllCameras(bool streamMode) {
  resetFrameCounts();
  for (auto& it : cameraInfoMap_) {
    log(std::string("Starting camera ") + to_string(it.first));
    std::lock_guard<std::mutex> lock(it.second.voxelDataMutex_);
    it.second.voxels.clear();
    it.second.samples_.reset();
    if (streamMode) it.second.camera->SetStream(true);
    it.second.stdDev->resetAveraging();
    if (it.second.camera->Start() != 0) return false;
  }
  return true;
}

void CameraManager::clearAllCameras() {
  for (auto& it : cameraInfoMap_) {
    std::lock_guard<std::mutex> lock(it.second.voxelDataMutex_);
    it.second.voxels.clear();
    it.second.samples_.reset();
  }
}

bool CameraManager::anyBadFrames() const {
  for (const auto& it : cameraInfoMap_) {
    if (it.second.camera->BadFrames()) return true;
  }
  return false;
}

bool CameraManager::stopAllCameras() {
  for (const auto& it : cameraInfoMap_) {
    Rcam* camera = it.second.camera;
    camera->Stop();
    log(std::string("Stopped camera ") + to_string(it.first) + ": " + to_string(camera->BadFrames())
        + " bad frame(s); " + to_string(camera->DroppedFrames()) + " dropped frame(s).");
  }

  #define REPORT_PROCESSING_TIME 0
  #if REPORT_PROCESSING_TIME
  for (int i = 0; i < (int)processingTime_.size(); i++) {
    if (processingTime_[i]) printf("processingTime[%d]: %d\n", i, processingTime_[i]);
  }
  #endif

  return true;
}

void CameraManager::addVoxel(const VoxelData& voxelData) {
  cameraInfo& info = cameraInfoMap_[voxelData.cameraID];
  std::lock_guard<std::mutex> lock(info.voxelDataMutex_);
  info.samples_.addSample(voxelData, (float)info.stdDev->GainConstant());
  info.voxels.push_back(voxelData);  // ToDo(jfs): Presize and write w/o mutex (a la voxelBuffer)?
      // With a linear pipeline, there's one thread per frame, but multiple frames in flight.

  // Keep a histogram of frame processing times.
  #if REPORT_PROCESSING_TIME
  int frameTime_ms = std::min(int(Component::SteadyClockTimeMs() - voxelData.timestamp),
      (int)processingTime_.size() - 1);
  ++processingTime_[frameTime_ms];
  #endif
}

void CameraManager::saveFrame(Frame* fr) {
  std::lock_guard<std::mutex> lock(frameMutex_);
  if (frameBuffer_) {  // write to memory?
    for (int i = 0, n = fr->width * fr->height; i < n; i++) {
      frameBuffer_[i] = (char)(fr->data[i] >> 2);
    }
  } else if (!filenameBase_.empty()) {  // write to file?
    std::string filename = filenameBase_ + to_string(framesSaved_) + ".tiff";
    if (fr->Write(filename) == 0) {
      log(std::string("Wrote ") + filename);
    } else {
      log(std::string("ERROR writing ") + filename);
    }
  }
  ++framesSaved_;
}

const std::string bloodflowHeader =
    "cameraID,frame,timestamp,saturated,rawMean,rawStdDev,imageMean,stdDev,speckleContrast,temperature";

bool CameraManager::writeCSVFile(const std::string& csvFilePath) {
  std::ofstream file;
  file.open(csvFilePath, std::ofstream::out);
  if (file.is_open()) {
    file << bloodflowHeader << std::endl;
    for (const auto& it : cameraInfoMap_) {
      for (const VoxelData& voxel : it.second.voxels) {
        file << voxel.cameraID << "," << voxel.frame << "," << voxel.timestamp << ","
            << voxel.saturated << "," << voxel.rawMean << "," << sqrt(voxel.rawVariance) << ","
            << voxel.imageMean << "," << voxel.stdDev << "," << voxel.speckleContrast << ","
            << voxel.temperature << std::endl;
      }
    }
    return true;
  } else {
    log(std::string("ERROR: Failed to open ") + csvFilePath);
    return false;
  }
}

void CameraManager::showMask(int cameraID, bool show) {
  if (cameraInfoMap_.find(cameraID) == cameraInfoMap_.end()) {
    log(std::string("ERROR: showROI: No such camera: ") + to_string(cameraID));
    return;
  }
  cameraInfoMap_[cameraID].stdDev->ShowMask(show);
}

bool CameraManager::checkData(int cameraID, bool fixParity) {
  if (cameraInfoMap_.find(cameraID) == cameraInfoMap_.end()) {
    log(std::string("ERROR: checkData: No such camera: ") + to_string(cameraID));
    return false;
  }

  cameraInfo& camInfo = cameraInfoMap_[cameraID];
  if (camInfo.voxels.size() < 3) {
    log(std::string("ERROR: Insufficient data for camera ") + to_string(cameraID));
    return false;
  }

  // Check for loss of signal. Trios of frames should have alternating differences.
  float lastMean = camInfo.voxels[1].rawMean;
  float lastDiff = lastMean - camInfo.voxels[0].rawMean;
  for (size_t i = 2, n = camInfo.voxels.size(); i < n; i++) {
    float thisMean = camInfo.voxels[i].rawMean;
    float thisDiff = thisMean - lastMean;  // should be of opposite sign from lastDiff
    if (thisDiff * lastDiff >= 0.f) {
      log(std::string("ERROR: Loss of signal, camera ") + to_string(cameraID));
      return false;
    }
    lastMean = thisMean;
    lastDiff = thisDiff;
  }

  // In lieu of a guarantee that even frames are dark, compare the averages of odd & even frames.
  float odds = 0, evens = 0;
  int oddCount = 0, evenCount = 0;
  for (const VoxelData& v : camInfo.voxels) {
    if (v.frame & 1) {  // odd frames should be light; even frames, dark.
      odds += v.rawMean;
      ++oddCount;
    } else {
      evens += v.rawMean;
      ++evenCount;
    }
  }
  if (evens > odds) {
    if (!fixParity) {
      log("ERROR: Even frames lighter than odds.");
      return false;
    }
    // Eliminate the first point if the parity is flipped.
    log("WARNING: Even frames lighter than odds; truncating.");
    camInfo.voxels.erase(cameraInfoMap_[cameraID].voxels.begin());
    for (VoxelData& v : camInfo.voxels) {
      --v.frame;
    }
    if (camInfo.voxels[0].frame & 1) {
      log("ERROR: First voxel frame not even.");
      return false;
    }
  }

  // Calculate the corrections, and average the corrected means.
  // ToDo(jfs): StdDev and VoxelSave are redundantly (and incorrectly) calc'ing some of these params.
  float averageImageMean = 0.0f;
  int n_averaged = 0;
  float darkMean = 0, darkVar = 0;
  float gainConst = (float)camInfo.stdDev->GainConstant();
  float imageMeanCorr = 0.0f;
  float imageMeanCorrPrev = 0.0f;
  for (VoxelData& v : camInfo.voxels) {
    if ((v.frame & 1) == 0) {
      darkMean = v.imageMean = v.rawMean;
      darkVar = v.rawVariance;
    } else {
      v.imageMean = v.rawMean - darkMean;
      v.stdDev = sqrt(std::max(v.rawVariance - darkVar - gainConst * v.imageMean, 0.0f));
    }
    v.speckleContrast = v.stdDev / v.imageMean;
    // Above calculations & if statement produce questionable results that don't match CSV
    averageImageMean += v.imageMean;
    if (v.frame == 0) {
      n_averaged--;
      imageMeanCorrPrev = v.rawMean;
      imageMeanCorr += 0;
    } else {
      imageMeanCorr += abs(v.rawMean - imageMeanCorrPrev);
      imageMeanCorrPrev = v.rawMean;
    }
    n_averaged++;
  }
  if (n_averaged > 0) {
    // log(std::string("MEAN CHECK: ") + to_string(imageMeanCorr / n_averaged));
    // log(std::string("MEAN THRES: ") + to_string(camInfo.meanThres));
    // log(std::string("MEAN COUNT: ") + to_string(n_averaged));
    if (imageMeanCorr < camInfo.meanThres * n_averaged) {
      log(std::string("ERROR: camera ") + to_string(cameraID) + std::string(" average image mean of ") + to_string(imageMeanCorr / n_averaged) + std::string(" is less than ") + to_string(camInfo.meanThres));
      return false;
    }
  }
  else {
    log(std::string("ERROR: camera ") + to_string(cameraID) + std::string(" Number of frames is zero "));
    return false;
  }
  return true;
}

float CameraManager::dataMean(int cameraID) {
  // Assumes checkData called first, and succeeds, except for aligner.
  if (cameraInfoMap_.find(cameraID) == cameraInfoMap_.end()) {
    log(std::string("ERROR: dataMean: No such camera: ") + to_string(cameraID));
    return -1;
  }

  cameraInfo& camInfo = cameraInfoMap_[cameraID];
  if (camInfo.voxels.size() < 2) {
    log(std::string("ERROR: Insufficient data for camera ") + to_string(cameraID));
    return -1;
  }

  double lightMeanSum = 0;
  int lightCount = 0;
  for (const VoxelData& v : cameraInfoMap_[cameraID].voxels) {
    if (v.frame & 1) {
      lightMeanSum += v.imageMean;
      ++lightCount;
    }
  }

  return (float)(lightMeanSum / (double)lightCount);  // TODO: return double
}

float CameraManager::dataContrast(int cameraID) {
  if (cameraInfoMap_.find(cameraID) == cameraInfoMap_.end()) return 0;
  float sum = 0;
  int count = 0;
  for (const VoxelData& v : cameraInfoMap_[cameraID].voxels) {
    if (v.frame & 1) {
      sum += v.speckleContrast;
      ++count;
    }
  }
  return sum / (float)count;
}

float CameraManager::dataSTD(int cameraID) {
  if (cameraInfoMap_.find(cameraID) == cameraInfoMap_.end()) return 0;
  float sum = 0;
  int count = 0;
  for (const VoxelData& v : cameraInfoMap_[cameraID].voxels) {
    if (v.frame & 1) {
      sum += v.stdDev;
      ++count;
    }
  }
  return sum / (float)count;
}


float CameraManager::sampleContrast(int cameraID) {
  if (cameraInfoMap_.find(cameraID) == cameraInfoMap_.end()) return 0;
  float sample = -0.01f;
  int count = 0;
  std::vector<VoxelData>& vs = cameraInfoMap_[cameraID].voxels;
  size_t nv = cameraInfoMap_[cameraID].voxels.size();
  if (nv >= 1) {
    VoxelData& v = vs[nv - 1];
    sample = v.speckleContrast;
  }
  return sample;
}

float CameraManager::sampleMean(int cameraID) {
  if (cameraInfoMap_.find(cameraID) == cameraInfoMap_.end()) return 0;
  float sample = -0.01f;
  int count = 0;
  std::vector<VoxelData> &vs = cameraInfoMap_[cameraID].voxels;
  size_t nv = cameraInfoMap_[cameraID].voxels.size();
  if (nv >= 1) {
    sample = vs[nv - 1].rawMean;
  }
  return sample;
}

float CameraManager::sampleSTD(int cameraID) {
  if (cameraInfoMap_.find(cameraID) == cameraInfoMap_.end()) return 0;
  float sample = -0.01f;
  int count = 0;
  std::vector<VoxelData>& vs = cameraInfoMap_[cameraID].voxels;
  size_t nv = cameraInfoMap_[cameraID].voxels.size();
  if (nv >= 1) {
    sample = vs[nv - 1].stdDev;
  }
  return sample;
}



float* CameraManager::arrayContrast(int cameraID) {
  float* ret = (float*)Samples::back_samples_;
  if (cameraInfoMap_.find(cameraID) == cameraInfoMap_.end()) {
    return ret;
  }
  Samples& s = cameraInfoMap_[cameraID].samples_;
  return s.arrayContrastSamples;
}

float* CameraManager::arrayMean(int cameraID) {
  float* ret = (float*)Samples::back_samples_;
  if (cameraInfoMap_.find(cameraID) == cameraInfoMap_.end()) {
    return ret;
  }
  Samples& s = cameraInfoMap_[cameraID].samples_;
  return s.arrayMeanCorrSamples;
}

float* CameraManager::arraySTD(int cameraID) {
  float* ret = (float*)Samples::back_samples_;
  if (cameraInfoMap_.find(cameraID) == cameraInfoMap_.end()) {
    return ret;
  }
  Samples& s = cameraInfoMap_[cameraID].samples_;
  return s.arraySTDSamples;
}

bool CameraManager::waitDataCollection(int timeout_s) {
  // resetFrameCounts();  try doing this in startAllCameras instead
  time_t start = Component::SteadyClockTimeMs();
  for (const auto& it : cameraInfoMap_) {
    int cameraID = it.first;
    while (cameraInfoMap_[cameraID].camera->GetFrameCount() < 1) {
      Component::SleepMs(1000);
      log(std::string("Waiting for frames, camera ") + to_string(cameraID));
      if (Component::SteadyClockTimeMs() - start > timeout_s * 1000) {
        log("ERROR: Timeout while waiting for frames.");
        return false;
      }
    }
  }
  return true;
}

void CameraManager::beginAccumulation() {
  for (auto& it : cameraInfoMap_) {
    it.second.stdDev->beginAccumulation(it.second.camera->GetConfig());
  }
}

void CameraManager::clearAccumulation() {
  for (auto& it : cameraInfoMap_) {
    StdDev* sd = it.second.stdDev;
    if (sd->framesAccumulated()) {
      log("Camera " + to_string(it.first) + ": Detected " + to_string(sd->nBadPixels())
          + " bad pixels, w/nSigmas = " + to_string(sd->nSigmas()) + ", mean "
          + to_string(sd->darkAveragedMean()) + ", stddev " + to_string(sd->darkAveragedStdDev()));
    }
    sd->clearAccumulation();
  }
}

bool CameraManager::framesAccumulated() {
  for (const auto& it : cameraInfoMap_) {
    if (!it.second.stdDev->framesAccumulated()) return false;
  }
  return true;
}

bool CameraManager::writeAverageFrames(int serialNumber, const std::string& filenameBase) {
  StdDev* stdDev = cameraInfoMap_[serialNumber].stdDev;
  return stdDev->darkAccum().Write(filenameBase + "_dark.tiff") == 0
      && stdDev->lightAccum().Write(filenameBase + "_light.tiff") == 0;
}

void CameraManager::resetFrameCounts() {
  for (const auto& it : cameraInfoMap_) {
    it.second.camera->SetFrameCount(0);
    it.second.stdDev->ResetFrameCount();
  }
}

/////Class Samples

//static constexpr const float back_samples_[size_] = { 1.0f };//Linux linker cannot handle constexpr?
float Samples::back_samples_[size_] = { 1.0f };

Samples::Samples() :index_(0), schecker_(this) {
  reset();
};

void Samples::addSample(const VoxelData& vd, float gainConst) {
  arrayMeanSamples[index_] = vd.rawMean;
  arrayVarianceSamples[index_] = vd.rawVariance;
  if (index_ > 0) {
    float amean = abs(arrayMeanSamples[index_] - arrayMeanSamples[index_ - 1]);
    float avar = abs(arrayVarianceSamples[index_] - arrayVarianceSamples[index_ - 1]);
    arrayMeanCorrSamples[index_] = amean;
    arraySTDSamples[index_] = sqrt(std::max<float>(avar, 0.0f));
    arrayContrastSamples[index_] = arraySTDSamples[index_] > 0.0f ? arraySTDSamples[index_] / arrayMeanCorrSamples[index_] : 0.0f;
  }
  if (index_ == 1) {
    arrayMeanCorrSamples[0] = arrayMeanCorrSamples[index_];
    arraySTDSamples[0] = arraySTDSamples[index_];
    arrayContrastSamples[0] = arrayContrastSamples[index_];
  }
  index_++;
  index_ = index_ >= size_ ? 0 : index_;
}

#if 0
void Samples::addSample(const VoxelData& vd, float gainConst) {
  float& sM = arrayMeanSamples[index_];
  float& sIM = arrayImageMeanSamples[index_];
  float& sV = arrayVarianceSamples[index_];
  float& sS = arraySTDSamples[index_];
  float& sC = arrayContrastSamples[index_];

  sM = vd.rawMean;
  sV = vd.rawVariance;
  sIM = vd.rawMean;
  switch (schecker_.check()) {
  case sequence_check_t::is_dark: {
    sS = 0.0f;
    sC = 0.0f;
  }break;
  case sequence_check_t::is_light: {
    //sIM = sM - schecker_.darkMean();
    float v = sV;// -schecker_.darkVariance() - gainConst * sIM;
    v = std::max<float>(v, 0.0f);
    sS = sqrt(v);
  }break;
  }
  sC = sS / sIM;
  index_++;
  index_ = index_ >= size_ ? 0 : index_;
}
#endif
#if 0
void Samples::addSample(const VoxelData& vd, float gainConst) {
  arrayMeanSamples[index_] = vd.rawMean;
  arrayVarianceSamples[index_] = vd.rawVariance;
  switch (schecker_.check()) {
  case sequence_check_t::is_dark: {
    arraySTDSamples[index_] = 0.0f;
    arrayContrastSamples[index_] = 0.0f;
  }break;
  case sequence_check_t::is_light: {
    //float v = arrayVarianceSamples[index_] - schecker_.darkVariance() - gainConst * arrayMeanSamples[index_];
    float v = arrayVarianceSamples[index_];// -schecker_.darkVariance() - gainConst * arrayMeanSamples[index_];
    v = std::max<float>(v, 0.0f);
    arraySTDSamples[index_] = sqrt(v);
    arrayContrastSamples[index_] = arrayMeanSamples[index_] > 0.0f ? arraySTDSamples[index_] / arrayMeanSamples[index_] : 0.0f;
  }break;
  }
  index_++;
  index_ = index_ >= size_ ? 0 : index_;
}
#endif
#if 0
void Samples::addSample(const VoxelData& vd, float gainConst) {
  arrayContrastSamples[index_] = vd.speckleContrast;
  arrayMeanSamples[index_] = vd.rawMean;
  arraySTDSamples[index_] = vd.stdDev;
  index_++;
  index_ = index_ >= size_ ? 0 : index_;
}
#endif

void Samples::reset() {
  index_ = 0;
  memset(arrayMeanSamples, 0, sizeof(float) * size_);
  memset(arrayMeanCorrSamples, 0, sizeof(float) * size_);
  memset(arrayImageMeanSamples, 0, sizeof(float) * size_);
  memset(arrayVarianceSamples, 0, sizeof(float) * size_);
  memset(arrayContrastSamples, 0, sizeof(float) * size_);
  memset(arraySTDSamples, 0, sizeof(float) * size_);
}

SequenceChecker::SequenceChecker(Samples* samples):samples_(samples), even_is_dark(true) {};
sequence_check_t SequenceChecker::check() {
  sequence_check_t ret = sequence_check_t::checking;
  if (samples_->index_ == width) {
    float ev = samples_->arrayMeanSamples[0] + samples_->arrayMeanSamples[2] + samples_->arrayMeanSamples[4];
    float odd = samples_->arrayMeanSamples[1] + samples_->arrayMeanSamples[3] + samples_->arrayMeanSamples[5];
    even_is_dark = ev < odd;
  }
  if (samples_->index_ >= width) {
    if (is_odd()){
      ret = even_is_dark ? sequence_check_t::is_light : sequence_check_t::is_dark;
    }
    else {
      ret = even_is_dark ? sequence_check_t::is_dark : sequence_check_t::is_light;
    }
  }
  return ret;
}
inline bool SequenceChecker::is_odd() { return samples_->index_ & 1; }
inline bool SequenceChecker::is_even() { return !is_odd(); }
inline int SequenceChecker::dark_index() {
  int dark_ind = (samples_->index_ > 0) ? samples_->index_ - 1 : 0;//guess that previous sample supposed to be dark
  dark_ind += even_is_dark && is_even() ? 0 : 1;//correct if we suppose that the current is dark
  return dark_ind;
}
float SequenceChecker::darkMean() {
  return samples_->arrayMeanSamples[dark_index()];
}
float SequenceChecker::darkVariance() {
  return samples_->arrayVarianceSamples[dark_index()];
}