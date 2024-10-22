#pragma once
#include "basic_structure_data.h"

class FullStructureData : public BasicStructureData {
 public:
  bool checkStructure() const override {
    if (false == BasicStructureData::checkStructure()) {
      return false;
    }

    //
    return true;
  }

  rbmd::Id* _h_molecules_id;
  /// bond
  rbmd::Id* _h_bond_type;
  rbmd::Id* _h_bond_id0;
  rbmd::Id* _h_bond_id1;
  
  //special
  rbmd::Real* _h_special_weights;
  rbmd::Id* _h_special_ids;
  rbmd::Id* _h_special_offsets;
  rbmd::Id* _h_special_offset_count;

  //
 rbmd::Id*  _h_atoms_vec_gro;
 rbmd::Id*  _h_countVector;

  /// angle
  rbmd::Id* _h_angle_type;
  rbmd::Id* _h_angle_id0;
  rbmd::Id* _h_angle_id1;
  rbmd::Id* _h_angle_id2;

  /// dihedral
  rbmd::Id* _h_dihedral_type;
  rbmd::Id* _h_dihedral_id0;
  rbmd::Id* _h_dihedral_id1;
  rbmd::Id* _h_dihedral_id2;
  rbmd::Id* _h_dihedral_id3;

  //rbmd::Id* _h_special_source_array;
  //rbmd::Id* _h_special_offsets_array;
};