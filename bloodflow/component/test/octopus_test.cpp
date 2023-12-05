#include "bloodflow/component/driver.h"
#include "bloodflow/component/octopus.h"

#include "googletest/googlemock/include/gmock/gmock.h"
#include "googletest/googletest/include/gtest/gtest.h"

using testing::Return;
using testing::_;

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

// Test subclass using MockFX3
class OctopusT : public Octopus {
 public:
  OctopusT(MockDriver* driver): Octopus(driver) {}
  ~OctopusT() { fx3_ = nullptr; }  // Don't delete MockDriver on stack.
};

class TestOctopus : public ::testing::Test {
 public:
  TestOctopus() {}
};

TEST(TestOctopus, OpenFailsIfDriverOpenFails) {
  MockDriver driver;
  OctopusT octo(&driver);
  EXPECT_CALL(driver, Open(0x4F12, 0)).Times(1).WillOnce(Return(-1));
  ASSERT_EQ(-1, octo.Open(0));
}

// ToDo(jfs): More!
