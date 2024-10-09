#pragma once
#include "object.h"

class System : public Object {
 public:
  virtual int Evolve() = 0;
};