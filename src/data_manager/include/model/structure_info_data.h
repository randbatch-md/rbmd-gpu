#pragma once

#include "../common/types.h"

class StructureInfoData {
 public:
  StructureInfoData() : _num_atoms1(0){};

 public:
  rbmd::Id* _num_atoms;
  rbmd::Id _num_atoms1;
  rbmd::Id* _num_bonds;
  rbmd::Id* _num_angles;
  rbmd::Id* _num_dihedrals;
  rbmd::Id* _num_impropers;
  rbmd::Id* _num_atoms_type;
  rbmd::Id* _num_bounds_type;
  rbmd::Id* _num_angles_type;
  rbmd::Id* _num_dihedrals_type;

  rbmd::Range* _range;
};
