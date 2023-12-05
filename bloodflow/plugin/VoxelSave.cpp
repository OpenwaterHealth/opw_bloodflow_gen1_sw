#include <string>
#include <iostream>
#include <algorithm>

#include "VoxelData.h"
#include "VoxelSave.h"

#include "bloodflow/component/frame.h"
#include "bloodflow/component/stddev.h"

void VoxelSave::saveFrames(int nFrames) {
  std::lock_guard<std::mutex> lock(mutex_);
  nFrames_ = nFrames;
}

void* VoxelSave::Exec(void* data) {
  Frame* fr = (Frame*)data;

  // Return a frame via callback, if requested.
  {
    std::lock_guard<std::mutex> lock(mutex_);
    if (nFrames_ > 0) {
      // Assume that if we're saving more than one frame, it includes both light and dark.
      if (nFrames_ > 1 || (fr->seq & 1)) {
        saveFrame_(fr);  // Writes under a mutex in caller.
        --nFrames_;
      }
    }
  }

  const StdDev::Tag* tag = StdDev::GetTag(fr);
  if (tag) {
    VoxelData voxelData;
    voxelData.cameraID = fr->serialNumber;
    voxelData.frame = fr->seq;
    voxelData.timestamp = fr->timestamp_ms_;
    voxelData.rawMean = (float)tag->rawMean;
    voxelData.rawVariance = (float)tag->rawVariance;
    voxelData.temperature = (float)fr->temperature;
    voxelData.saturated = tag->saturated;
    //calculate the rest of the voxelData values, like in the CameraManager::checkData
    float darkMean = 0.0f, darkVar = 0.0f;
    if ((voxelData.frame & 1) == 0) {
      darkMean = voxelData.imageMean = voxelData.rawMean;
      darkVar = voxelData.rawVariance;
    }
    else {
      voxelData.imageMean = voxelData.rawMean - darkMean;
      voxelData.stdDev = sqrt(std::max(voxelData.rawVariance - darkVar - (float)gainConstant_ * voxelData.imageMean, 0.0f));
    }
    voxelData.speckleContrast = voxelData.stdDev / voxelData.imageMean; //!!! There will be NaN for 0 imageMean
    addVoxel_(voxelData);
  }//!!! handle null tag
  return data;
}
