#include "bloodflow/component/driver.h"
#include "bloodflow/component/rcam.h"
#include "bloodflow/plugin/CameraManager.h"
#include "googletest/googlemock/include/gmock/gmock.h"
#include "googletest/googletest/include/gtest/gtest.h"

using json = nlohmann::json;
using testing::Return;
using testing::_;

// Test metadata
const char* const TestJSON = R"foo({
  "cameraParameters": {
    "exposureTime_ms": 9.9,
    "resolutionY0": 0,
    "resolutionY1": 1024,
    "blackLevelCompensation": false,
    "gain": 1.0
  },
  "delayParameters": {
    "pulseWidth_s": 0.000100
  }
}
)foo";

class MockRcam: public Rcam {
 public:
  MockRcam(): Rcam(nullptr) {
    param_.bits = 16;
    param_.fpms = 1;
    param_.height = 2080;
    param_.width = 2712;
    memcpy(param_.model,"MOCK",4);
    param_.px_clk = 100000000;
  }
  MOCK_METHOD1(Open, int(int nth));
  MOCK_METHOD0(SerialNumber, int());
  MOCK_METHOD0(Reset, int());
  MOCK_METHOD0(Start, int());
  MOCK_METHOD0(Stop, void());
  MOCK_METHOD0(GetFrameCount, int());
  MOCK_METHOD1(SetExposure, void(double exp));
  MOCK_METHOD1(SetBLC, void(bool blc_on));
  MOCK_METHOD1(SetGain, void(double gain));
  MOCK_METHOD1(resize, void(size_t n));
  MOCK_METHOD1(SetStream, void(bool stream));
  MOCK_METHOD4(SubWindow2Point, void(int x0, int y0, int x1, int y1));
};

int NCameras = 0;
int Driver::NumDevices(uint16_t pid, bool recheck) { return NCameras; }

class TestCameraManager : public ::testing::Test {
 public:
  TestCameraManager() {}
};

class CameraManagerT : public CameraManager {
 public:
  static MockRcam* cam_;
  static Rcam* newRcam() { Rcam* cam = cam_; cam_ = nullptr; return cam; }  // hands off ownership
  CameraManagerT(): CameraManager(newRcam) { cam_ = new MockRcam(); }
  void setData(int cameraID, int n, float* v) {
    for(int i = 0; i < n; i++) {
      VoxelData voxel;
      voxel.frame = i;
      voxel.cameraID = cameraID;
      voxel.rawMean = v[i];
      addVoxel(voxel);
    }
  }
};
MockRcam* CameraManagerT::cam_;

TEST(TestCameraManager, CreateAndDestroyWork) {
  CameraManagerT cmgr;
  ASSERT_EQ(0, 0);
}

TEST(TestCameraManager, InitWorksWithNoCameras) {
  json systemParameters = json::parse(std::stringstream(TestJSON));
  CameraManagerT cmgr;
  ASSERT_TRUE(cmgr.init(systemParameters));
}

TEST(TestCameraManager, InitFailsIfCameraOpenFails) {
  json systemParameters = json::parse(std::stringstream(TestJSON));
  CameraManagerT cmgr;
  NCameras = 1;
  EXPECT_CALL((*cmgr.cam_), Open(0)).Times(1).WillOnce(Return(-1));
  ASSERT_FALSE(cmgr.init(systemParameters));
}

static void expectCameraTest(MockRcam& cam) {
  EXPECT_CALL(cam, Open(0)).Times(1).WillOnce(Return(0));
  EXPECT_CALL(cam, SerialNumber()).Times(2).WillRepeatedly(Return(101));
  EXPECT_CALL(cam, GetFrameCount()).Times(2).WillOnce(Return(0)).WillOnce(Return(4));
  EXPECT_CALL(cam, SetStream(true)).Times(1);  // stream on
  EXPECT_CALL(cam, Start()).Times(1).WillOnce(Return(0));  // stream on
  EXPECT_CALL(cam, Stop()).Times(1);  // stream on
  EXPECT_CALL(cam, SetStream(false)).Times(1);  // stream off
  EXPECT_CALL(cam, SetGain(1.)).Times(1);
  EXPECT_CALL(cam, SetBLC(false)).Times(1);
  EXPECT_CALL(cam, SetExposure(_)).Times(1);
  EXPECT_CALL(cam, resize(_)).Times(2);
  EXPECT_CALL(cam, SubWindow2Point(0, 0, 2712, 1024)).Times(1);
}

TEST(TestCameraManager, ReturnsCorrectSerialNumber) {
  json systemParameters = json::parse(std::stringstream(TestJSON));
  CameraManagerT cmgr;
  NCameras = 1;
  Rcam* savedRcam = cmgr.cam_;
  expectCameraTest(*cmgr.cam_);
  ASSERT_TRUE(cmgr.init(systemParameters));  // nulls out cmgr.cam_
  ASSERT_EQ(101, cmgr.serialNumber(0));
  ASSERT_EQ(savedRcam, cmgr.getCamera(101));
  ASSERT_EQ(0., cmgr.meanMean(101));
  ASSERT_EQ(0., cmgr.darkMean(101));
  ASSERT_EQ(0., cmgr.lightMean(101));
}

TEST(TestCameraManager, CheckDataFailsWithBadCameraID) {
  json systemParameters = json::parse(std::stringstream(TestJSON));
  CameraManagerT cmgr;
  expectCameraTest(*cmgr.cam_);
  ASSERT_TRUE(cmgr.init(systemParameters));
  ASSERT_FALSE(cmgr.checkData(1 /* no such camera*/, false));
}

TEST(TestCameraManager, CheckDataFailsWithNoData) {
  json systemParameters = json::parse(std::stringstream(TestJSON));
  CameraManagerT cmgr;
  NCameras = 1;
  expectCameraTest(*cmgr.cam_);
  ASSERT_TRUE(cmgr.init(systemParameters));  // nulls out cmgr.cam_
  ASSERT_FALSE(cmgr.checkData(101, false));
}

TEST(TestCameraManager, CheckDataWorks) {
  json systemParameters = json::parse(std::stringstream(TestJSON));
  CameraManagerT cmgr;
  NCameras = 1;
  expectCameraTest(*cmgr.cam_);
  ASSERT_TRUE(cmgr.init(systemParameters));
  float data[] = {0.f, 1.f, 3.f, 2.f};  // reverses direction
  cmgr.setData(101, 4, data);
  ASSERT_FALSE(cmgr.checkData(101, false));  // fails with bad data
  cmgr.clearAllCameras();
  float data2[] = {0.f, 1.f, 0.1f, 1.1f};  // same directions
  cmgr.setData(101, 4, data2);
  ASSERT_TRUE(cmgr.checkData(101, false));  // succeeds with good data
  cmgr.clearAllCameras();
  float data3[] = {1.f, 0.f, 1.f, 0.f};  // parity wrong
  cmgr.setData(101, 4, data3);
  ASSERT_FALSE(cmgr.checkData(101, false));  // fails w/fixParity false
  ASSERT_TRUE(cmgr.checkData(101, true));  // succeeds w/fixParity true
}

// ToDo(jfs): More ...
