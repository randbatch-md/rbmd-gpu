#pragma once
#include "../common/object.h"
// #include "system.h"
#include "ensemble.h"
#include "../common/types.h"
#include "json/reader.h"
#include "json/value.h"
#include "../output/include/output.h"
extern int test_current_step;

class Simulate : public Object {
 public:
  Simulate(std::shared_ptr<Ensemble>&simulate_pipeline, std::shared_ptr<Output>& output);
  virtual ~Simulate() = default;

 public:
  void Init();
  int Execute();
  bool KeepGoing();

 protected:
  // Json::Value _exec_node;
  std::shared_ptr<Ensemble>& _simulate_pipeline;
  std::shared_ptr<Output>& _output;
  float _time_step;
  float _current_time;

  int _num_steps;
  int _current_step;
};