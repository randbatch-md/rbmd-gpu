#pragma once

#include "../common/types.h"

class StructureInfoData {
public:
  ~StructureInfoData();

  rbmd::Id* _num_atoms = nullptr;
  rbmd::Id* _num_bonds = nullptr;
  rbmd::Id* _num_angles = nullptr;
  rbmd::Id* _num_dihedrals = nullptr;
  rbmd::Id* _num_impropers = nullptr;
  rbmd::Id* _num_atoms_type = nullptr;
  rbmd::Id* _num_bounds_type = nullptr;
  rbmd::Id* _num_angles_type = nullptr;
  rbmd::Id* _num_dihedrals_type = nullptr;
  rbmd::Id* _num_impropers_type = nullptr;
};

inline StructureInfoData::~StructureInfoData() {
  FREE(_num_atoms);
  FREE(_num_bonds);
  FREE(_num_angles);
  FREE(_num_dihedrals);
  FREE(_num_impropers);
  FREE(_num_atoms_type);
  FREE(_num_bounds_type);
  FREE(_num_angles_type);
  FREE(_num_dihedrals_type);
  FREE(_num_impropers_type);
}