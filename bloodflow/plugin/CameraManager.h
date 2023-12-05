#pragma once

#include <fstream>
#include <map>
#include <mutex>
#include <string>
#include <vector>
#include <algorithm>

#include "bloodflow/component/stddev.h"
#include "bloodflow/third_party/json-develop/single_include/nlohmann/json.hpp"

#include "VoxelData.h"

class Rcam;
class VoxelSave;

// Utility class
class SequenceChecker;
class Samples;
enum class sequence_check_t {
  checking,
  is_dark,
  is_light
};
class SequenceChecker {
  friend class Samples;
public:
  static const int width = 6;
public:
  SequenceChecker(Samples* samples);
  sequence_check_t check();
  inline bool is_odd();
  inline bool is_even();
  inline int dark_index();
  float darkMean();
  float darkVariance();
private:
  Samples* samples_;
  bool even_is_dark;
};
//A utility class containing arrays of current values available for real-time plotting
class Samples {
  friend class SequenceChecker;
public:
  static const int size_ = 1000;
  //static constexpr const float back_samples_[size_] = { 1.0f };//Linux linker cannot handle
  static float back_samples_[size_];
  Samples();

  void addSample(const VoxelData& vd, float gainConst);
  void reset();

  float arrayMeanSamples[size_];
  float arrayMeanCorrSamples[size_];
  float arrayImageMeanSamples[size_];
  float arrayVarianceSamples[size_];
  float arrayContrastSamples[size_];
  float arraySTDSamples[size_];
private:
  uint32_t index_;
  SequenceChecker schecker_;
};

// Encapsulate all the expertise on [multi]camera handling for scanning.
class CameraManager {
 public:
  CameraManager(Rcam*(*newRcam)(), void (*log)(const char*) = nullptr):
      newRcam_(newRcam), log_(log) {}
  virtual ~CameraManager();

  // Initializes [multi]cameras and processing pipelines for speckle contrast experiments
  virtual bool init(const nlohmann::json& systemParameters);

  // Test a camera.
  bool testCamera(Rcam*);

  // Reinit the ROIs, which has to happen before a 2nd calibration (or more).
  void ReinitROIs();

  // Return the serial number for the nth camera.
  int serialNumber(int nth);

  // Get a camera by serialNumber.
  Rcam* getCamera(int serialNumber) { return cameraInfoMap_[serialNumber].camera; }

  // Start & stop rcam objects
  bool startAllCameras(bool streamMode);
  bool stopAllCameras();

  // Return whether any cameras currently have a nonzero bad frame count.
  bool anyBadFrames() const;

  // Check the efficacy of the acquired data, optionally fixing parity if necessary.
  bool checkData(int cameraID, bool fixParity);

  // Write voxel buffer to csv file
  bool writeCSVFile(const std::string& csvFilePath);

  // Show the ROI (or not) for a given camera.
  void showMask(int serialNumber, bool show);

  // Clear voxels from all cameras.
  void clearAllCameras();

  // Return the mean of the voxel imageMeans & speckle contrasts, for a given camera.
  float dataMean(int cameraID);
  float dataContrast(int cameraID);
  float dataSTD(int cameraID);

  // Return current samples
  float sampleContrast(int cameraID);
  float sampleMean(int cameraID);
  float sampleSTD(int cameraID);
  float *arrayContrast(int cameraID);
  float *arrayMean(int cameraID);
  float *arraySTD(int cameraID);

  // Start image saving, for calibration.
  void saveImages(int serialNumber, int nFrames, char* buffer);

  // Wait for frames to be saved.
  bool waitForFrames(int nFrames);

  // Return the running means (per-camera).
  double meanMean(int serialNumber) { return cameraInfoMap_[serialNumber].stdDev->meanMean(); }
  double darkMean(int serialNumber) { return cameraInfoMap_[serialNumber].stdDev->darkMean(); }
  double lightMean(int serialNumber) { return cameraInfoMap_[serialNumber].stdDev->lightMean(); }

  // Frame accumulation calls, passed to all cameras
  void beginAccumulation(), clearAccumulation();
  bool framesAccumulated();

  // Write the averaged frames to filenameBase + {_dark,_light}.tiff.
  bool writeAverageFrames(int serialNumber, const std::string& filenameBase);

  // Wait for data collection to begin after button press. Return false on timeout.
  bool waitDataCollection(int timeout_s);

  // Zero all camera frame counts.
  void resetFrameCounts();

  // Send log messages here.
  void logWith(void (*log)(const char*)) { log_ = log; }

  // Log, via log_ if possible, stdout otherwise.
  void log(const char* s) const;
  void log(const std::string& s) const { log(s.c_str()); }

 protected:
  // Add processed data to voxel buffer (called from voxelSave thread)
  void addVoxel(const VoxelData& voxelData);

  // Save a frame (with mutex).
  void saveFrame(Frame* frame);

  // Rcam factory
  Rcam* (*newRcam_)();

  // Container mapping each cameraID attached [key: (int) cameraID#, value: struct]
  struct cameraInfo {
    int cameraID = -1;
    double meanThres = 5.0;
    Rcam* camera;
    StdDev* stdDev;
    VoxelSave* voxelSave;
    std::vector<VoxelData> voxels;
    Samples samples_;
    std::mutex voxelDataMutex_;
  };
  std::map<int, cameraInfo> cameraInfoMap_;

  // Histogram of frame processing times, in ms.
  std::vector<int> processingTime_ = std::vector<int>(1000);

  int framesSaved_ = 0;          // count of saveFrame() calls
  char* frameBuffer_ = nullptr;  // for saving frames to memory; copy frames to this buffer
  std::string filenameBase_;     // for saving frames to files (only one camera at a time, currently)
  std::mutex frameMutex_;        // libtiff is not reentrant
  void (*log_)(const char* s);
};
