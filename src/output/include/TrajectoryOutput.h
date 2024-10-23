#pragma once
#include "../executioner/executioner.h"
#include "output.h"
#include <iostream>
#include <fstream>

class TrajectoryOutput : public Output {
 public:
  TrajectoryOutput();
  virtual ~TrajectoryOutput() = default;

  void Init() override;
  void Execute() override;

  bool ShouldOutput();

 private:
  rbmd::Id _num_atoms;


  rbmd::Id _interval;
  rbmd::Id _atom_num;
};