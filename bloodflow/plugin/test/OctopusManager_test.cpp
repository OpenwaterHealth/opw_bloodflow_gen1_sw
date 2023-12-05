#include "bloodflow/plugin/OctopusManager.h"
#include "googletest/googlemock/include/gmock/gmock.h"
#include "googletest/googletest/include/gtest/gtest.h"

using json = nlohmann::json;
using testing::Return;
using testing::_;

// Octopus methods are not virtual yet, so bypass creating Mocktopus. Use a real Octopus and pass
// it a MockDriver.
class MockDriver : public Driver {
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

// Test subclass using MockDriver
class OctopusT : public Octopus {
 public:
  OctopusT(MockDriver* driver): Octopus(driver) {}
  ~OctopusT() { fx3_ = nullptr; }  // Don't delete MockDriver on stack.
};

class TestOctopusManager : public ::testing::Test {
 public:
  TestOctopusManager() {}
};

// Test metadata
const char* const TestJSON = R"foo({
  "AOMParameters": { "AOM4Volt_V": 0.26, "AOM4Freq_MHz": 100 },
  "cameraParameters": { "numImages": 2500, "frameAcquisitionRate_Hz": 28 },
  "delayParameters": {
    "pulseWidth_s": 0.0001,
    "triggerOffset_s": 0.01,
    "frameValidDelay_s": 0.009000000000000001,
    "frameValidWidth_s": 0.0025,
    "chEDelay_s": 0.0089
  },
  "hardwareParameters": { "hwTrigger": false }
}
)foo";

TEST(TestOctopusManager, CreateAndDestroyWork) {
  OctopusManager omgr(new Octopus(new MockDriver()));
  ASSERT_EQ(0, 0);
}

TEST(TestOctopusManager, InitWorks) {
  json systemParameters = json::parse(std::stringstream(TestJSON));
  MockDriver driver;
  EXPECT_CALL(driver, Open(_, _)).Times(1).WillOnce(Return(0));
  EXPECT_CALL(driver, Flash(_, _)).Times(1).WillOnce(Return(0));
  EXPECT_CALL(driver, SerialNumber()).Times(1).WillOnce(Return(42));
  ON_CALL(driver, DataIn(2, _)).WillByDefault(Return(2));
  ON_CALL(driver, DataOut(2, _)).WillByDefault(Return(2));
  ON_CALL(driver, DataOut(4, _)).WillByDefault(Return(4));
  OctopusManager omgr(new OctopusT(&driver));
  ASSERT_TRUE(omgr.init(systemParameters));
}

TEST(TestOctopusManager, InitWithCapturePeriodTooLongFails) {
  json systemParameters = json::parse(std::stringstream(TestJSON));
  systemParameters["cameraParameters"]["numImages"] = 5000;
  MockDriver driver;
  ON_CALL(driver, DataIn(2, _)).WillByDefault(Return(2));
  ON_CALL(driver, DataOut(2, _)).WillByDefault(Return(2));
  ON_CALL(driver, DataOut(4, _)).WillByDefault(Return(4));
  OctopusManager omgr(new OctopusT(&driver));
  ASSERT_FALSE(omgr.init(systemParameters));
}
