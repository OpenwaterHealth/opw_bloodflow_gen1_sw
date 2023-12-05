#include "bloodflow/component/execnode.h"
#include "bloodflow/component/time.h"
#include "googletest/googlemock/include/gmock/gmock.h"
#include "googletest/googletest/include/gtest/gtest.h"

using testing::Return;
using testing::_;

class MockNode : public ExecNode {
 public:
  MOCK_METHOD2(Schedule, bool(void* data, ExecNode* producer));
  MOCK_METHOD1(Execute, void(void* data));
  MOCK_METHOD1(resize, void(size_t n));
  MOCK_METHOD2(CleanUp, void(void* data, ExecNode* consumer));
};

class ExecNodeT : public ExecNode {  // for access to protected members
 public:
  ExecNodeT() {}
  ~ExecNodeT() {}
  std::vector<ExecNode*>& producers() { return producers_; }
  std::vector<ExecNode*>& consumers() { return consumers_; }
  CircularBuffer<std::pair<ExecNode*, void*>>& scheduledNodes() {
    return thread_manager_.scheduledNodes();
  }
  CircularBuffer<void*>& downQueue() { return down_queue_; }
  CircularBuffer<void*>& upQueue() { return up_queue_; }
  bool Schedule2(void* data, ExecNode* producer) { return Schedule(data, producer); }
  void Execute2(void* data) { Execute(data); }
  void CleanUp2(void* data, ExecNode* consumer) { CleanUp(data, consumer); }
  void Run2(const std::vector<bool>& run_list, void* data) { Run(run_list, data); }
  void* Exec(void* data) override { Component::SleepMs(100); return data; }  // don't run too fast!
};

class TestExecNode : public ::testing::Test {
 public:
  TestExecNode() {}
};

TEST(TestExecNode, CreateAndDestroyWork) {
  ExecNode node;
  ASSERT_EQ(0, 0);  // i.e., doesn't crash
}

TEST(TestExecNode, AddProducerWorks) {
  ExecNodeT e1, e2;
  e2.AddProducer(&e1);
  ASSERT_EQ(e1.consumers()[0], &e2);
  ASSERT_EQ(e2.producers()[0], &e1);
  ASSERT_TRUE(e1.IsRoot());
  ASSERT_FALSE(e1.IsLeaf());
  ASSERT_FALSE(e2.IsRoot());
  ASSERT_TRUE(e2.IsLeaf());
}

TEST(TestExecNode, ProduceCallsNodeSchedule) {
  ExecNode::SetNumberThreads(1);
  ExecNodeT e1;
  MockNode e2;
  e1.consumers().push_back(&e2);
  void* data = (void*)0xDEADCAF;
  EXPECT_CALL(e2, Schedule(data, &e1)).Times(1).WillOnce(Return(false));
  e1.Produce(data);
  ASSERT_TRUE(true);
}

#ifdef _MSC_VER  // this is flaky on linux (ToDo(jfs): Fix this.)
TEST(TestExecNode, ProduceCallsThreadManagerSchedule) {
  ExecNode::SetNumberThreads(1);
  ExecNodeT e1, e2;
  e2.AddProducer(&e1);
  void* data = (void*)0xDEADCAF;
  e1.Produce(data);
  ASSERT_EQ(e1.scheduledNodes().Peek(0).first, &e2);
  while (!e2.IsExecDone()) Component::SleepMs(100);
}
#endif

TEST(TestExecNode, ResizeCallsConsumer) {
  ExecNode::SetNumberThreads(1);
  ExecNodeT e1;
  MockNode e2;
  e2.AddProducer(&e1);
  EXPECT_CALL(e2, resize(30)).Times(1);
  e1.resize(30);
  ASSERT_EQ(30, e1.downQueue().size());
  ASSERT_EQ(30, e1.upQueue().size());
}

TEST(TestExecNode, ScheduleReturnsFalseForMultipleInputs) {
  ExecNodeT e1, e2, e3;
  e3.AddProducer(&e1);
  e3.AddProducer(&e2);
  ASSERT_FALSE(e3.Schedule2(nullptr, &e1));
}

TEST(TestExecNode, ExecuteLeafCallsProducerCleanup) {
  MockNode e1;
  ExecNodeT e2;
  e2.AddProducer(&e1);
  void* data = (void*)0xDEADCAF;
  e2.Schedule2(data, &e1);
  EXPECT_CALL(e1, CleanUp(data, &e2)).Times(1);
  e2.Execute2(data);
}

TEST(TestExecNode, ExecuteNonleafCallsConsumerSchedule) {
  ExecNodeT e1;
  MockNode e2;
  e2.AddProducer(&e1);
  void* data = (void*)0xDEADCAF;
  e1.Schedule2(data, nullptr);
  EXPECT_CALL(e2, Schedule(data, &e1)).Times(1);
  e1.Execute2(data);
}

TEST(TestExecNode, CleanUpCallsProducerCleanUp) {
  MockNode e1;
  ExecNodeT e2;
  e2.AddProducer(&e1);
  void* data = (void*)0xDEADCAF;
  e2.Schedule2(data, &e1);
  EXPECT_CALL(e1, CleanUp(data, &e2)).Times(1);
  e2.CleanUp2(data, nullptr);
}

TEST(TestExecNode, RunWithSingleConsumerCallsExecute) {
  ExecNodeT e1;
  MockNode e2;
  e2.AddProducer(&e1);
  void* data = (void*)0xDEADCAF;
  EXPECT_CALL(e2, Execute(data)).Times(1);
  std::vector<bool> run_list;
  run_list.push_back(true);
  e1.Run2(run_list, data);
}
