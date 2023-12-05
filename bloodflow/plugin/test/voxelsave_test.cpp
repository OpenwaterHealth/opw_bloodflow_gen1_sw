#include "bloodflow/component/stddev.h"
#include "bloodflow/plugin/VoxelData.h"
#include "bloodflow/plugin/VoxelSave.h"
#include "googletest/googlemock/include/gmock/gmock.h"
#include "googletest/googletest/include/gtest/gtest.h"

using testing::Return;
using testing::_;

class TestVoxelSave : public ::testing::Test {
 public:
  TestVoxelSave() {}
};

class VoxelSaveT : public VoxelSave {
 public:
  VoxelSaveT(std::function<void(const VoxelData&)> addVoxel): VoxelSave(addVoxel, nullptr, 0.0) {}
  void Exec2(void* data) { Exec(data); }
};

TEST(TestVoxelSave, ExecWorks) {
  int addVoxelCalls = 0;
  VoxelData voxelData;
  VoxelSaveT voxelSave([&](const VoxelData& v){ ++addVoxelCalls; voxelData = v; });
  Frame frame;
  StdDev::Tag tag;
  tag.rawMean = 8.;
  tag.rawVariance = 16.;
  frame.AddTag(&tag);
  voxelSave.Exec2(&frame);
  ASSERT_EQ(1, addVoxelCalls);
  ASSERT_EQ(8, voxelData.rawMean);
  ASSERT_EQ(16.f, voxelData.rawVariance);
}
