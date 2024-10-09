#pragma once
#include <vector>

#include "../common/object.h"
#include "../common/types.h"

class ForceFieldData : public Object {
 public:
  virtual bool checkForceField() const = 0;
};
