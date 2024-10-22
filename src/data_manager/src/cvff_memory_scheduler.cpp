#include "./include/scheduler/cvff_memory_scheduler.h"

bool CVFFMemoryScheduler::asyncMemoryH2D() {
  if (false == MemoryScheduler::asyncMemoryH2D()) {
    // log
    return false;
  }

  auto& num_atoms_type = *(_structure_info_data->_num_atoms_type);
  auto& num_atoms = *(_structure_info_data->_num_atoms);
  auto& num_bonds = *(_structure_info_data->_num_bonds);
  auto& num_bonds_type = *(_structure_info_data->_num_bounds_type);
  auto& num_angles = *(_structure_info_data->_num_angles);
  auto& num_angles_type = *(_structure_info_data->_num_angles_type);
  auto& num_dihedrals = *(_structure_info_data->_num_dihedrals);
  auto& num_dihedrals_type = *(_structure_info_data->_num_dihedrals_type);
  auto sd = std::dynamic_pointer_cast<FullStructureData>(_structure_data);
  auto fd = std::dynamic_pointer_cast<CVFFForceFieldData>(_force_field_data);

  /// copy data
  _device_data->_d_molecular_id.resize(num_atoms);
  _device_data->_d_bond_type.resize(num_bonds);
  _device_data->_d_bond_id0.resize(num_bonds);
  _device_data->_d_bond_id1.resize(num_bonds);
  _device_data->_d_angle_type.resize(num_angles);
  _device_data->_d_angle_id0.resize(num_angles);
  _device_data->_d_angle_id1.resize(num_angles);
  _device_data->_d_angle_id2.resize(num_angles);
  _device_data->_d_dihedral_type.resize(num_dihedrals);
  _device_data->_d_dihedral_id0.resize(num_dihedrals);
  _device_data->_d_dihedral_id1.resize(num_dihedrals);
  _device_data->_d_dihedral_id2.resize(num_dihedrals);
  _device_data->_d_dihedral_id3.resize(num_dihedrals);
  _device_data->_d_charge.resize(num_atoms);
  _device_data->_d_special_source_array.resize(*(sd->_num_special_source_array));
  _device_data->_d_special_offsets_array.resize(*(sd->_num_special_offsets_array));
  _device_data->_d_special_weights.resize(*(sd->_num_special_weights));
  _device_data->_d_special_ids.resize(*(sd->_num_special_ids));
  _device_data->_d_special_offsets.resize(*(sd->_num_special_offsets));
  _device_data->_d_special_count.resize(sizeof(sd->_h_special_offset_count));

  thrust::copy(sd->_h_special_source_array, sd->_h_special_source_array + *(sd->_num_special_source_array),
      _device_data->_d_special_source_array.begin());
  thrust::copy(sd->_h_special_offsets_array, sd->_h_special_offsets_array + *(sd->_num_special_offsets_array),
      _device_data->_d_special_offsets_array.begin());

  /// charge
  thrust::copy(sd->_h_charge, sd->_h_charge + num_atoms,
      _device_data->_d_charge.begin());

  /// molecular id
  thrust::copy(sd->_h_molecules_id, sd->_h_molecules_id + num_atoms,
               _device_data->_d_molecular_id.begin());
  /// bond
  thrust::copy(sd->_h_bond_type, sd->_h_bond_type + num_bonds,
               _device_data->_d_bond_type.begin());
  thrust::copy(sd->_h_bond_id0, sd->_h_bond_id0 + num_bonds,
               _device_data->_d_bond_id0.begin());
  thrust::copy(sd->_h_bond_id1, sd->_h_bond_id1 + num_bonds,
               _device_data->_d_bond_id1.begin());
  //special
  thrust::copy(sd->_h_special_weights, sd->_h_special_weights + *(sd->_num_special_weights),
      _device_data->_d_special_weights.begin());
  thrust::copy(sd->_h_special_ids, sd->_h_special_ids + *(sd->_num_special_ids),
      _device_data->_d_special_ids.begin());
  thrust::copy(sd->_h_special_offsets, sd->_h_special_offsets + *(sd->_num_special_offsets),
      _device_data->_d_special_offsets.begin());
  thrust::copy(sd->_h_special_offset_count, sd->_h_special_offset_count+
    sizeof(sd->_h_special_offset_count),_device_data->_d_special_count.begin());

  // atoms_vec  sourceArray
  thrust::copy(sd->_h_atoms_vec_gro, sd->_h_atoms_vec_gro + sizeof(sd->_h_atoms_vec_gro),
    _device_data->_d_atoms_vec.begin());

  thrust::copy(sd->_h_countVector, sd->_h_countVector + sizeof(sd->_h_countVector),
    _device_data->_d_atoms_count.begin());

  //atoms_vec offsetArray
  //cpu上运行
  std::vector<rbmd::Id> cumulative_offsets;
  cumulative_offsets.push_back(0); // 初始偏移量为0

  for (int i = 0; i < sizeof(sd->_h_countVector); ++i)
  {
    cumulative_offsets.push_back(cumulative_offsets.back() + sd->_h_countVector[i]);
  }
  thrust::copy(cumulative_offsets.begin(), cumulative_offsets.end(),
  _device_data->_d_atoms_offset.begin());

  //GPU上运行
  // thrust::device_vector<rbmd::Id> d_atoms_offset_temp(sd->_h_countVector.size());
  // thrust::copy(sd->_h_countVector.begin(), sd->_h_countVector.end(),
  //   d_atoms_offset_temp.begin());
  //
  // _device_data->_d_atoms_offset.resize(d_atoms_offset_temp.size() + 1);
  // _device_data->_d_atoms_offset[0] = 0;
  // thrust::exclusive_scan(d_atoms_offset_temp.begin(), d_atoms_offset_temp.end(),
  //   _device_data->_d_atoms_offset.begin() + 1);


  /// angle
  thrust::copy(sd->_h_angle_type, sd->_h_angle_type + num_angles,
               _device_data->_d_angle_type.begin());
  thrust::copy(sd->_h_angle_id0, sd->_h_angle_id0 + num_angles,
               _device_data->_d_angle_id0.begin());
  thrust::copy(sd->_h_angle_id1, sd->_h_angle_id1 + num_angles,
               _device_data->_d_angle_id1.begin());
  thrust::copy(sd->_h_angle_id2, sd->_h_angle_id2 + num_angles,
               _device_data->_d_angle_id2.begin());

  /// dihedral
  thrust::copy(sd->_h_dihedral_type, sd->_h_dihedral_type + num_dihedrals,
               _device_data->_d_dihedral_type.begin());
  thrust::copy(sd->_h_dihedral_id0, sd->_h_dihedral_id0 + num_dihedrals,
               _device_data->_d_dihedral_id0.begin());
  thrust::copy(sd->_h_dihedral_id1, sd->_h_dihedral_id1 + num_dihedrals,
               _device_data->_d_dihedral_id1.begin());
  thrust::copy(sd->_h_dihedral_id2, sd->_h_dihedral_id2 + num_dihedrals,
               _device_data->_d_dihedral_id2.begin());
  thrust::copy(sd->_h_dihedral_id3, sd->_h_dihedral_id3 + num_dihedrals,
               _device_data->_d_dihedral_id3.begin());

  /// copy force field
  _device_data->_d_eps.resize(num_atoms_type);
  _device_data->_d_mass.resize(num_atoms_type);
  _device_data->_d_sigma.resize(num_atoms_type);
  _device_data->_d_bond_coeffs_k.resize(num_bonds_type);
  _device_data->_d_bond_coeffs_equilibrium.resize(num_bonds_type);
  _device_data->_d_angle_coeffs_k.resize(num_angles_type);
  _device_data->_d_angle_coeffs_equilibrium.resize(num_angles_type);
  _device_data->_d_dihedral_coeffs_k.resize(num_dihedrals_type);
  _device_data->_d_dihedral_coeffs_sign.resize(num_dihedrals_type);
  _device_data->_d_dihedral_coeffs_multiplicity.resize(num_dihedrals_type);

  /// eps
  thrust::copy(fd->_h_eps, fd->_h_eps + num_atoms_type,
               _device_data->_d_eps.begin());

  /// mass
  thrust::copy(fd->_h_mass, fd->_h_mass + num_atoms_type,
               _device_data->_d_mass.begin());

  /// sigma
  thrust::copy(fd->_h_sigma, fd->_h_sigma + num_atoms_type,
               _device_data->_d_sigma.begin());

  /// bond
  thrust::copy(fd->_h_bond_coeffs_k, fd->_h_bond_coeffs_k + num_bonds_type,
               _device_data->_d_bond_coeffs_k.begin());
  thrust::copy(fd->_h_bond_coeffs_equilibrium,
               fd->_h_bond_coeffs_equilibrium + num_bonds_type,
               _device_data->_d_bond_coeffs_equilibrium.begin());

  /// angle
  thrust::copy(fd->_h_angle_coeffs_k, fd->_h_angle_coeffs_k + num_angles_type,
               _device_data->_d_angle_coeffs_k.begin());
  thrust::copy(fd->_h_angle_coeffs_equilibrium,
               fd->_h_angle_coeffs_equilibrium + num_angles_type,
               _device_data->_d_angle_coeffs_equilibrium.begin());

  /// dihedral
  thrust::copy(fd->_h_dihedral_coeffs_k,
               fd->_h_dihedral_coeffs_k + num_dihedrals_type,
               _device_data->_d_dihedral_coeffs_k.begin());
  thrust::copy(fd->_h_dihedral_coeffs_sign,
               fd->_h_dihedral_coeffs_sign + num_dihedrals_type,
               _device_data->_d_dihedral_coeffs_sign.begin());
  thrust::copy(fd->_h_dihedral_coeffs_multiplicity,
               fd->_h_dihedral_coeffs_multiplicity + num_dihedrals_type,
               _device_data->_d_dihedral_coeffs_multiplicity.begin());
  return true;
}

bool CVFFMemoryScheduler::asyncMemoryD2H() { return true; }
