#pragma once
#include "ensemble.h"

class NVEensemble : public Ensemble {
 public:
  NVEensemble();
  virtual ~NVEensemble() = default;

 protected:
  void Init() override;
  void Presolve() override;
  void Solve() override;
  void Postsolve() override;

 private:
};