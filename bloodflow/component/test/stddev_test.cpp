#include "bloodflow/component/stddev.h"
#include "googletest/googlemock/include/gmock/gmock.h"
#include "googletest/googletest/include/gtest/gtest.h"

using testing::Return;
using testing::_;

class TestStddev : public ::testing::Test {
 public:
  TestStddev() {}
};

class StdDevT : public StdDev {
 public:
  StdDevT() {}
  void Exec2(void* data) { Exec(data); }
  void AtExit2(void* data) { AtExit(data); }
};

TEST(TestStddev, ExecWorks) {
  StdDevT stddev;
  Frame frame;
  uint16_t pixels[] = { 1, 2, 3, 4 };
  frame.data = pixels;
  frame.bits = 16;
  frame.width = 2;
  frame.height = 2;
  stddev.resize(1);
  stddev.Exec2(&frame);
  ASSERT_EQ(2.5, StdDev::GetTag(&frame)->rawMean);
  ASSERT_DOUBLE_EQ(1.1180339887498949, sqrt(StdDev::GetTag(&frame)->rawVariance));
  ASSERT_EQ(2.5, stddev.meanMean());
  ASSERT_EQ(2.5, stddev.darkMean());  // ToDo(jfs): Feed this thing more frames to average.
  ASSERT_EQ(0, stddev.lightMean());
}

TEST(TestStddev, BeginAndEndAccumulationWork) {
  StdDevT stddev;
  Frame& darkAccum = stddev.darkAccum();
  Frame& lightAccum = stddev.lightAccum();
  ASSERT_EQ(darkAccum.width, 0);
  ASSERT_EQ(darkAccum.height, 0);
  ASSERT_EQ(lightAccum.width, 0);
  ASSERT_EQ(lightAccum.height, 0);
  ASSERT_FALSE(stddev.framesAccumulated());

  Frame frame(2, 3);
  stddev.beginAccumulation(&frame);
  ASSERT_EQ(darkAccum.width, 2);
  ASSERT_EQ(darkAccum.height, 3);
  ASSERT_EQ(lightAccum.width, 2);
  ASSERT_EQ(lightAccum.height, 3);
  for (int i = 0; i < 6; i++) {
    ASSERT_EQ(darkAccum.data[i], 0);
    ASSERT_EQ(lightAccum.data[i], 0);
  }
  ASSERT_FALSE(stddev.framesAccumulated());  // ToDo(jfs): accumulate some frames.

  stddev.clearAccumulation();
  ASSERT_EQ(darkAccum.width, 0);
  ASSERT_EQ(darkAccum.height, 0);
  ASSERT_EQ(lightAccum.width, 0);
  ASSERT_EQ(lightAccum.height, 0);
  ASSERT_FALSE(stddev.framesAccumulated());
}

TEST(TestStddev, AveragingWorks) {
  StdDevT stddev;
  Frame frame;
  uint16_t pixels[] = { 1, 2, 3, 4 };
  frame.data = pixels;
  frame.bits = 16;
  frame.width = 2;
  frame.height = 2;
  frame.seq = 0;
  stddev.resize(1);
  stddev.Exec2(&frame);
  stddev.AtExit2(&frame);
  ASSERT_EQ(2.5, stddev.meanMean());
  pixels[0] = 0;  // dark frame
  stddev.Exec2(&frame);
  stddev.AtExit2(&frame);
  ASSERT_EQ(2.375, stddev.darkMean());
  pixels[0] = 3;
  frame.seq = 1;
  stddev.Exec2(&frame);
  stddev.AtExit2(&frame);
  ASSERT_EQ(3.0, stddev.lightMean());
  stddev.resetAveraging();
  ASSERT_EQ(0., stddev.darkMean());
  ASSERT_EQ(0., stddev.lightMean());
}

TEST(TestStddev, ROIWorks) {
  StdDevT stddev;
  Frame frame;
  uint16_t pixels[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };  // 3x3 image
  frame.data = pixels;
  frame.bits = 10;
  frame.width = 3;
  frame.height = 3;
  stddev.resize(3);
  stddev.Exec2(&frame);
  ASSERT_EQ(5, StdDev::GetTag(&frame)->rawMean);
  ASSERT_DOUBLE_EQ(2.5819888974716112, sqrt(StdDev::GetTag(&frame)->rawVariance));
  frame.ClearTags();
  pixels[0] = pixels[2] = pixels[6] = pixels[8] = 1023;
  stddev.Exec2(&frame);
  ASSERT_DOUBLE_EQ(457.44444444444446, StdDev::GetTag(&frame)->rawMean);
  ASSERT_EQ(4, StdDev::GetTag(&frame)->saturated);
  frame.ClearTags();
  stddev.SetROI(1, 1, 1, 3, 3);  // exclude corner pixels
  stddev.Exec2(&frame);
  ASSERT_EQ(5, StdDev::GetTag(&frame)->rawMean);
  ASSERT_DOUBLE_EQ(2, sqrt(StdDev::GetTag(&frame)->rawVariance));
}
