#include <thread>

#include "bloodflow/component/frame.h"
#include "bloodflow/component/TiffInterface.h"
#include "bloodflow/third_party/tiff/inc/tiffio.h"
#include "bloodflow/third_party/tiff/inc/tiff.h"
#include "googletest/googlemock/include/gmock/gmock.h"
#include "googletest/googletest/include/gtest/gtest.h"

using testing::Return;
using testing::_;

// Mock a TiffInterface for testing Frame's tiff i/o (see gMock docs for details).
class MockTiff : public TiffInterface {
 public:
  MOCK_METHOD2(Open, bool(const char*, const char*));
  MOCK_METHOD0(Close, void());
  MOCK_METHOD2(GetField, int(ttag_t tag, int* vp));
  MOCK_METHOD2(GetField, int(ttag_t tag, uint16_t* vp));
  MOCK_METHOD2(SetField, int(ttag_t tag, int value));
  MOCK_METHOD2(ReadScanline, int(tdata_t buf, uint32 row));
  MOCK_METHOD2(WriteScanline, int(tdata_t buf, uint32 row));
  MOCK_METHOD0(WriteDirectory, int());
};

// Subclass Frame using MockTiff.
class FrameTest : public Frame {
 public:
  FrameTest(MockTiff* tiff): Frame(20, 10) { tiff_ = tiff; /* happens after Init() */ }
  FrameTest() {}
  ~FrameTest() { tiff_ = NULL; }  // keep parent class from deleting stack object
  bool hasOutput() { return tiff_ != NULL; }
  TiffInterface* getTiff() { return tiff_; }
  void DoInitTiff() { InitTiff(); }
};

class TestFrame : public ::testing::Test {
 public:
  TestFrame() {}
};

// Read TEST(x, name) as "test that <name>".
// And a line such as
//   EXPECT_CALL(mockTiff, Open(_, _)).Times(1).WillOnce(Return(false));
// means, expect a call to mockTiff.Open(), with any args; when it comes, return false (i.e., err).

TEST(TestFrame, FrameBitsIs16) {
  // Note the default is 16, but when Rcam sends frames they'll be 10-bit.
  Frame f1(3, 2);
  ASSERT_EQ(16, f1.bits);
  Frame f2;
  ASSERT_EQ(0, f2.bits);
}

TEST(TestFrame, FrameAssignsCorrectly) {
  Frame f1(20, 10);
  f1[2][3] = 23, f1[3][2] = 32;
  Frame f2(f1);  // test copy constructor
  ASSERT_EQ(f2.width, f1.width);
  ASSERT_EQ(f2.height, f1.height);
  ASSERT_EQ(f2[2][3], 23);
  ASSERT_EQ(f2[3][2], 32);
  Frame f3;  // N.B.: Frame f3 = f1 does NOT do the same as below.
  f3 = f1;  // test copy by assignment
  ASSERT_EQ(f3.width, f1.width);
  ASSERT_EQ(f3.height, f1.height);
  ASSERT_EQ(f3[2][3], 23);
  ASSERT_EQ(f3[3][2], 32);
}

TEST(TestFrame, FrameCreatesDefaultOutput) {
  FrameTest test;
  test.DoInitTiff();  // should create tiff_ TiffInterface
  ASSERT_TRUE(test.hasOutput());
}

TEST(TestFrame, FramesCreateSeparateOutputs) {  // Necessary for tiff i/o to be reentrant.
  FrameTest test1;
  FrameTest test2;
  test1.DoInitTiff();  // These call InitTiff, and should create two TiffInterfaces.
  test2.DoInitTiff();  // (Previously, tiff_ was static.)
  ASSERT_FALSE(test1.getTiff() == test2.getTiff());
}

TEST(TestFrame, FrameDoesNotCopyOutput) {
  testing::NiceMock<MockTiff> mockTiff;
  FrameTest test1(&mockTiff);
  FrameTest test2(test1);  // this should not copy tiff_
  ASSERT_TRUE(test1.getTiff() == &mockTiff);
  ASSERT_TRUE(test2.getTiff() == NULL);
}

TEST(TestFrame, FramesGetDistinctTimestamps) {
  Frame f1, f2;
  f1.SetTimestamp();  // stamps f1 with 'now'
  std::this_thread::sleep_for(std::chrono::milliseconds(10));
  f2.SetTimestamp();  // stamps f2 with 'now'
  ASSERT_TRUE(f2.timestamp_ms_ > f1.timestamp_ms_);
}

TEST(TestFrame, WriteFailsWhenOpenFails) {
  testing::NiceMock<MockTiff> mockTiff;
  FrameTest test(&mockTiff);
  EXPECT_CALL(mockTiff, Open(_, _)).Times(1).WillOnce(Return(false));
  ASSERT_EQ(-1, test.Write("foo"));
}

#if 0  // not worth it; WriteScanline will fail if disk is full.
TEST(TestFrame, WriteFailsWhenSetFieldFails) {
  testing::NiceMock<MockTiff> mockTiff;
  FrameTest test(&mockTiff);
  EXPECT_CALL(mockTiff, Open(_, _)).Times(1).WillOnce(Return(true));
  EXPECT_CALL(mockTiff, SetField(_, _)).Times(1).WillOnce(Return(-1));
  ASSERT_EQ(-1, test.Write("foo"));
}
#endif

TEST(TestFrame, WriteFailsWhenWriteScanlineFails) {
  testing::NiceMock<MockTiff> mockTiff;
  FrameTest test(&mockTiff);
  EXPECT_CALL(mockTiff, Open(_, _)).Times(3).WillRepeatedly(Return(true));  // now let Open pass
  EXPECT_CALL(mockTiff, SetField(_, _)).Times(36).WillRepeatedly(Return(1));  // these are ok
  EXPECT_CALL(mockTiff, WriteScanline(_, 0)).Times(3).WillRepeatedly(Return(-1));  // this fails
  ASSERT_EQ(-1, test.Write("foo"));
  test.bits = 10;
  ASSERT_EQ(-1, test.Write("foo"));
  test.bits = 8;
  ASSERT_EQ(-1, test.Write("foo"));
}

TEST(TestFrame, WriteFailsWhenWriteDirectoryFails) {
  testing::NiceMock<MockTiff> mockTiff;
  FrameTest test(&mockTiff);
  EXPECT_CALL(mockTiff, Open(_, _)).Times(1).WillOnce(Return(true));
  EXPECT_CALL(mockTiff, SetField(_, _)).Times(12).WillRepeatedly(Return(1));
  EXPECT_CALL(mockTiff, WriteScanline(_, _)).Times(10).WillRepeatedly(Return(1));
  EXPECT_CALL(mockTiff, WriteDirectory()).Times(1).WillOnce(Return(-1));  // now this one fails
  ASSERT_EQ(-1, test.Write("foo"));
}

TEST(TestFrame, RealTiffWrites) {
  Frame f(32, 24);  // no subclass; really Write()s to files
  ASSERT_EQ(0, f.Write("foo16"));  // default 16-bit
  f.bits = 8;  // Try saving in 8-bit format also.
  ASSERT_EQ(0, f.Write("foo8"));
  f.bits = 10;  // Rcam default 10-bit
  ASSERT_EQ(0, f.Write("foo10"));
}

TEST(TestFrame, AddAndGetTag) {
  struct TestTag : Frame::Tag {
    int i;
  } tag;
  Frame fr;
  fr.AddTag(&tag);  // add a tag and make sure we can get it back
  ASSERT_EQ(fr.GetTag<TestTag>(), &tag);
}

TEST(TestFrame, TagNullWhenNoTag) {
  struct TestTag : Frame::Tag {
    int i;
  };
  Frame fr;
  ASSERT_EQ(fr.GetTag<TestTag>(), nullptr);
}

TEST(TestFrame, GetCorrectTag) {
  struct TestTag1 : Frame::Tag {
    int i;
  } tag1;
  struct TestTag2 : Frame::Tag {
    double i;
  } tag2;
  Frame fr;
  fr.AddTag(&tag1);  // add tags of two types and make sure we get them back right
  EXPECT_EQ(fr.GetTag<TestTag1>(), &tag1);
  EXPECT_EQ(fr.GetTag<TestTag2>(), nullptr);
  fr.AddTag(&tag2);
  EXPECT_EQ(fr.GetTag<TestTag2>(), &tag2);
}

TEST(TestFrame, ClearTagRemovesTags) {
  struct TestTag : Frame::Tag {
    int i;
  } tag;
  Frame fr;
  fr.AddTag(&tag);
  EXPECT_EQ(fr.GetTag<TestTag>(), &tag);
  fr.ClearTags();  // this should clear the tag
  EXPECT_EQ(fr.GetTag<TestTag>(), nullptr);
}
