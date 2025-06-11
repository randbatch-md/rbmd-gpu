#pragma once

#include "common/types.h"
#include "force_field_data.h"

class MACEForceFieldData : public ForceFieldData {
 public:
  bool checkForceField() const override { return true; }

 public:
  /// mass
  //rbmd::Real* _h_mass;

};
