#include "bloodflow/component/time.h"

#include "googletest/googlemock/include/gmock/gmock.h"
#include "googletest/googletest/include/gtest/gtest.h"

class TestTime : public ::testing::Test {
 public:
  TestTime() {}
};

TEST(TestTime, SteadyClockTimeMsIncrementsAfterSleep) {
  time_t t1 = Component::SteadyClockTimeMs();
  Component::SleepMs(10);
  time_t t2 = Component::SteadyClockTimeMs();
  ASSERT_GE(t2 - t1, 10);
}
