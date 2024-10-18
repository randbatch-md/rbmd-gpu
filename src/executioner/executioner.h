#pragma once
#include "common/object.h"
// #include "system.h"
#include "../simulate_pipeline/include/ensemble.h"
#include "common/types.h"
#include "json/reader.h"
#include "json/value.h"

extern int test_current_step;

class Executioner : public Object {
 public:
  Executioner(/*const Json::Value& node, */ std::shared_ptr<Ensemble>&
                  simulate_pipeline);
  virtual ~Executioner() = default;

 public:
  void Init();
  int Execute();
  bool KeepGoing();

 protected:
  // Json::Value _exec_node;
  std::shared_ptr<Ensemble>& _simulate_pipeline;

  float _time_step;
  float _current_time;

  int _num_steps;
  int _current_step;
};