#pragma once

#include <functional>
#include <mutex>

#include "bloodflow/component/execnode.h"

class Frame;
struct VoxelData;

// Node to fill in voxels from earlier calculations and forward to owner via callback.
class VoxelSave : public ExecNode {
 public:
  VoxelSave(std::function<void(const VoxelData&)> addVoxel, std::function<void(Frame*)> saveFrame, double gainConstant)
      : addVoxel_(addVoxel), saveFrame_(saveFrame), gainConstant_(gainConstant){}
  ~VoxelSave() {};

  // Set up for file saving, via the saveFrame callback.
  void saveFrames(int nFrames);

 protected:
  void* Exec(void* data) override;
  std::function<void(const VoxelData&)> addVoxel_;
  std::function<void(Frame*)> saveFrame_;

  int nFrames_ = 0;  // number of frames to save
  std::mutex mutex_;

  double gainConstant_;
};
