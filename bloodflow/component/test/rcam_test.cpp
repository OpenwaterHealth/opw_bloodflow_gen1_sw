#include "bloodflow/component/driver.h"
#include "bloodflow/component/rcam.h"
#include "googletest/googlemock/include/gmock/gmock.h"
#include "googletest/googletest/include/gtest/gtest.h"

using testing::Return;
using testing::_;

class MockDriver: public Driver {
 public:
  MOCK_METHOD2(Open, int(uint16_t pid, int n));
  MOCK_METHOD0(Reset, int());
  MOCK_METHOD0(SerialNumber, int());
  MOCK_METHOD0(Close, void());
  MOCK_METHOD2(Flash, int(int len, uint8_t* data));
  MOCK_METHOD5(CmdWrite,
      int(uint8_t req, uint16_t wValue, uint16_t wIndex, uint16_t len, uint8_t* buf));
  MOCK_METHOD5(CmdRead,
      int(uint8_t req, uint16_t wValue, uint16_t wIndex, uint16_t len, uint8_t* buf));
  MOCK_METHOD2(DataIn, int(int len, uint8_t* data));
  MOCK_METHOD2(DataOut, int(int len, uint8_t* data));
  MOCK_METHOD3(BulkInBuffers, int(int buflen, int n_buffers, int timeout_ms));
  MOCK_METHOD0(BulkInStart, int());
  MOCK_METHOD2(BulkInData, int(int& len, uint8_t*& data));
  MOCK_METHOD0(BulkInStop, void());
  MOCK_METHOD0(Reattach, void());
};
int Driver::NumDevices(uint16_t pid, bool recheck) { return 0; }

// Test subclass using MockFX3
class RcamTest: public Rcam {
 public:
  RcamTest(MockDriver* driver): Rcam(driver) {}
  ~RcamTest() { fx3_ = nullptr; }  // don't delete MockDriver on stack
  const Frame& FrCfg() { return fr_cfg_; }
  RcamParam& param() { return param_; }
};

class TestRcam: public ::testing::Test {
 public:
  TestRcam() {}
};

TEST(TestRcam, OpenFailsIfFX3OpenFails) {
  testing::NiceMock<MockDriver> mockDriver;
  RcamTest test(&mockDriver);
  EXPECT_CALL(mockDriver, Open(Rcam::pid(), 0)).Times(1).WillOnce(Return(-1));
  ASSERT_EQ(-1, test.Open(0));
}

TEST(TestRcam, OpenFailsIfFlashFails) {
  testing::NiceMock<MockDriver> mockDriver;
  RcamTest test(&mockDriver);
  EXPECT_CALL(mockDriver, Open(Rcam::pid(), 0)).Times(1).WillOnce(Return(0));
  EXPECT_CALL(mockDriver, Flash(_, _)).Times(1).WillOnce(Return(-1));
  ASSERT_EQ(-1, test.Open(0));
}

TEST(TestRcam, OpenFailsIfParamReadFails) {
  testing::NiceMock<MockDriver> mockDriver;
  RcamTest test(&mockDriver);
  EXPECT_CALL(mockDriver, Open(Rcam::pid(), 0)).Times(1).WillOnce(Return(0));
  EXPECT_CALL(mockDriver, Flash(_, _)).Times(1).WillOnce(Return(0));
  EXPECT_CALL(mockDriver, CmdRead(_, _, _, _, _)).Times(3).WillOnce(Return(0));
  ASSERT_EQ(-1, test.Open(0));
}

TEST(TestRcam, OpenFailsIfNoCameraName) {
  testing::NiceMock<MockDriver> mockDriver;
  RcamTest test(&mockDriver);
  EXPECT_CALL(mockDriver, Open(Rcam::pid(), 0)).Times(1).WillOnce(Return(0));
  EXPECT_CALL(mockDriver, Flash(_, _)).Times(1).WillOnce(Return(0));
  EXPECT_CALL(mockDriver, CmdRead(_, _, _, _, _)).Times(1).WillOnce(Return(sizeof(RcamParam)));
  ASSERT_EQ(-1, test.Open(0));
}

static int FakeCmdRead(uint8_t req, uint16_t wValue, uint16_t wIndex, uint16_t len, uint8_t* buf) {
  RcamParam* param = (RcamParam*)buf;
#if defined(_MSC_VER)
  strcpy_s(param->model, RCAM_PARAM_MODEL_TEXT_ZISE, "Hi, Max!");
#else
  //This GCC compiler does not have secure version of strcpy_s 
  strncpy(param->model, "Hi, Max!", RCAM_PARAM_MODEL_TEXT_ZISE);
#endif
  return sizeof(RcamParam);
}

TEST(TestRcam, OpenFailsIfBulkInBuffersFails) {
  testing::NiceMock<MockDriver> mockDriver;
  RcamTest test(&mockDriver);
  EXPECT_CALL(mockDriver, Open(Rcam::pid(), 0)).Times(1).WillOnce(Return(0));
  EXPECT_CALL(mockDriver, Flash(_, _)).Times(1).WillOnce(Return(0));
  ON_CALL(mockDriver, CmdRead(_, _, _, _, _)).WillByDefault(FakeCmdRead);
  EXPECT_CALL(mockDriver, BulkInBuffers(_, _, _)).Times(1).WillOnce(Return(-1));
  ASSERT_EQ(-1, test.Open(0));
}

TEST(TestRcam, SubWindow2PointSetsCfg) {
  testing::NiceMock<MockDriver> mockDriver;
  RcamTest test(&mockDriver);
  test.param().width = 1024;
  test.param().height = 1024;
  EXPECT_CALL(mockDriver, CmdWrite(_, _, _, _, _)).Times(6);  // just make sure this gets called
  test.SubWindow2Point(2, 3, 999, 474);
  ASSERT_EQ(test.FrCfg().width, 997);
  ASSERT_EQ(test.FrCfg().height, 471);
}

TEST(TestRcam, ResetFailsIfDriverResetFails) {
  testing::NiceMock<MockDriver> mockDriver;
  RcamTest test(&mockDriver);
  EXPECT_CALL(mockDriver, Reset()).Times(1).WillOnce(Return(-1));
  ASSERT_EQ(-1, test.Reset());
}

TEST(TestRcam, ResetFailsIfDriverFlashFails) {
  testing::NiceMock<MockDriver> mockDriver;
  RcamTest test(&mockDriver);
  EXPECT_CALL(mockDriver, Reset()).Times(1).WillOnce(Return(0));
  EXPECT_CALL(mockDriver, Flash(_, _)).Times(1).WillOnce(Return(-1));
  ASSERT_EQ(-1, test.Reset());
}

TEST(TestRcam, ResetSucceeds) {
  testing::NiceMock<MockDriver> mockDriver;
  RcamTest test(&mockDriver);
  test.param().width = 1024;
  test.param().height = 1024;
  test.SubWindow2Point(0, 0, 1024, 1024);
  EXPECT_CALL(mockDriver, Reset()).Times(1).WillOnce(Return(0));
  EXPECT_CALL(mockDriver, Flash(_, _)).Times(1).WillOnce(Return(0));
  EXPECT_CALL(mockDriver, CmdWrite(_, _, _, _, _)).Times(7).WillRepeatedly(Return(4));
  ASSERT_EQ(0, test.Reset());
}

// TODO(jfs): Hire consultant; more here ...
