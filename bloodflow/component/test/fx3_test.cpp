// fx3 unit test
// Note: w/o a virtual cypress-fx3 interface, it's hard to write a full test. CyFX3 is slated to
// go away, so this is judged not worth it. Nor is the windows container on the cloud going to have
// ow-fx3-driver.

#include "bloodflow/component/fx3.h"

#include "googletest/googlemock/include/gmock/gmock.h"
#include "googletest/googletest/include/gtest/gtest.h"

class TestFX3 : public ::testing::Test {
 public:
  TestFX3() {}
};

TEST(TestFX3, NumDevicesIsZero) {
  ASSERT_EQ(0, FX3::NumDevices(0x4F10));  // happens to be the Rcam pid
}

TEST(TestFX3, OpenFails) {
  FX3 fx3;
  ASSERT_EQ(-1, fx3.Open(0x4F10, 0));
}

// Nothing else works due to assertions.
