#pragma once

// Voxel data to write to CSV
struct VoxelData {
  int cameraID = 0;
  int frame = 0;
  int saturated = 0;
  time_t timestamp = 0;
  float rawMean = 0;
  float rawVariance = 0;
  float imageMean = 0;
  float stdDev = 0;
  float speckleContrast = 0;
  float temperature = 0;
};
