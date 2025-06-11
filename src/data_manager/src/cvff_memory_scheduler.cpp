#include "./include/scheduler/cvff_memory_scheduler.h"
#include "data_manager.h"

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
  auto& num_impropers = *(_structure_info_data->_num_impropers);
  auto& num_impropers_type = *(_structure_info_data->_num_impropers_type);
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
  _device_data->_d_angle_id_vec.resize(num_angles);

  _device_data->_d_dihedral_type.resize(num_dihedrals);
  _device_data->_d_dihedral_id0.resize(num_dihedrals);
  _device_data->_d_dihedral_id1.resize(num_dihedrals);
  _device_data->_d_dihedral_id2.resize(num_dihedrals);
  _device_data->_d_dihedral_id3.resize(num_dihedrals);

  _device_data->_d_improper_type.resize(num_impropers);
  _device_data->_d_improper_id0.resize(num_impropers);
  _device_data->_d_improper_id1.resize(num_impropers);
  _device_data->_d_improper_id2.resize(num_impropers);
  _device_data->_d_improper_id3.resize(num_impropers);

  _device_data->_d_charge.resize(num_atoms);

  _device_data->_d_atoms_vec.resize(sd->_num_atoms_vec_gro);
  _device_data->_d_atoms_count.resize(sd->_num_count_vector);
  _device_data->_d_atoms_offset.resize(sd->_num_atoms_offset );

  _device_data->_d_special_weights.resize(sd->_num_special_weights);
  _device_data->_d_special_ids.resize(sd->_num_special_ids);
  _device_data->_d_special_offsets.resize(sd->_num_special_offsets);
  _device_data->_d_special_count.resize(sd->_num_special_offset_count);

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

  //special weights and ids
  thrust::copy(sd->_h_special_weights, sd->_h_special_weights + sd->_num_special_weights,
      _device_data->_d_special_weights.begin());
  thrust::copy(sd->_h_special_ids, sd->_h_special_ids + sd->_num_special_ids,
      _device_data->_d_special_ids.begin());
  thrust::copy(sd->_h_special_offsets, sd->_h_special_offsets + sd->_num_special_offsets,
      _device_data->_d_special_offsets.begin());
  thrust::copy(sd->_h_special_offset_count, sd->_h_special_offset_count +
    sd->_num_special_offset_count,_device_data->_d_special_count.begin());

  // special atoms_vec
  thrust::copy(sd->_h_atoms_vec_gro, sd->_h_atoms_vec_gro + sd->_num_atoms_vec_gro,
      _device_data->_d_atoms_vec.begin());

  thrust::copy(sd->_h_count_vector, sd->_h_count_vector + sd->_num_count_vector,
      _device_data->_d_atoms_count.begin());

  thrust::copy(sd->_h_atoms_offset, sd->_h_atoms_offset + sd->_num_atoms_offset,
  _device_data->_d_atoms_offset.begin());

  //GPU
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
  thrust::copy(sd->_h_angle_id_vec, sd->_h_angle_id_vec + num_angles,
               _device_data->_d_angle_id_vec.begin());
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

  //improper
  thrust::copy(sd->_h_improper_type, sd->_h_improper_type + num_impropers,
             _device_data->_d_improper_type.begin());
  thrust::copy(sd->_h_improper_id0, sd->_h_improper_id0 + num_impropers,
               _device_data->_d_improper_id0.begin());
  thrust::copy(sd->_h_improper_id1, sd->_h_improper_id1 + num_impropers,
               _device_data->_d_improper_id1.begin());
  thrust::copy(sd->_h_improper_id2, sd->_h_improper_id2 + num_impropers,
               _device_data->_d_improper_id2.begin());
  thrust::copy(sd->_h_improper_id3, sd->_h_improper_id3 + num_impropers,
               _device_data->_d_improper_id3.begin());

  /// (2) copy force field
  /// mass
  // _device_data->_d_mass.resize(num_atoms_type);
  // thrust::copy(fd->_h_mass, fd->_h_mass + num_atoms_type,
  //              _device_data->_d_mass.begin());
  /// eps
 _device_data->_d_eps.resize(num_atoms_type);
  thrust::copy(fd->_h_eps, fd->_h_eps + num_atoms_type,
               _device_data->_d_eps.begin());
  /// sigma
  _device_data->_d_sigma.resize(num_atoms_type);
  thrust::copy(fd->_h_sigma, fd->_h_sigma + num_atoms_type,
               _device_data->_d_sigma.begin());
  /// bond
  _device_data->_d_bond_coeffs_k.resize(num_bonds_type);
  _device_data->_d_bond_coeffs_equilibrium.resize(num_bonds_type);
  thrust::copy(fd->_h_bond_coeffs_k, fd->_h_bond_coeffs_k + num_bonds_type,
               _device_data->_d_bond_coeffs_k.begin());
  thrust::copy(fd->_h_bond_coeffs_equilibrium,
               fd->_h_bond_coeffs_equilibrium + num_bonds_type,
               _device_data->_d_bond_coeffs_equilibrium.begin());
  /// angle
  _device_data->_d_angle_coeffs_k.resize(num_angles_type);
  _device_data->_d_angle_coeffs_equilibrium.resize(num_angles_type);
  thrust::copy(fd->_h_angle_coeffs_k, fd->_h_angle_coeffs_k + num_angles_type,
               _device_data->_d_angle_coeffs_k.begin());
  thrust::copy(fd->_h_angle_coeffs_equilibrium,
               fd->_h_angle_coeffs_equilibrium + num_angles_type,
               _device_data->_d_angle_coeffs_equilibrium.begin());
  //dihedral
  std::string dihedral_type = "null";

  if(*(_structure_info_data->_num_dihedrals)) {
    dihedral_type = DataManager::getInstance().getConfigData()->
      Get<std::string>("dihedral_type", "hyper_parameters", "force_field");
  }
  if (dihedral_type == "harmonic") {
    _device_data->_d_dihedral_coeffs_k.resize(num_dihedrals_type);
    _device_data->_d_dihedral_coeffs_sign.resize(num_dihedrals_type);
    _device_data->_d_dihedral_coeffs_multiplicity.resize(num_dihedrals_type);

    thrust::copy(fd->_h_dihedral_coeffs_k,
           fd->_h_dihedral_coeffs_k + num_dihedrals_type,
           _device_data->_d_dihedral_coeffs_k.begin());
    thrust::copy(fd->_h_dihedral_coeffs_sign,
                 fd->_h_dihedral_coeffs_sign + num_dihedrals_type,
                 _device_data->_d_dihedral_coeffs_sign.begin());
    thrust::copy(fd->_h_dihedral_coeffs_multiplicity,
                 fd->_h_dihedral_coeffs_multiplicity + num_dihedrals_type,
                 _device_data->_d_dihedral_coeffs_multiplicity.begin());
  }
  else if (dihedral_type == "opls") {
    _device_data->_d_dihedral_coeffs_k1.resize(num_dihedrals_type);
    _device_data->_d_dihedral_coeffs_k2.resize(num_dihedrals_type);
    _device_data->_d_dihedral_coeffs_k3.resize(num_dihedrals_type);
    _device_data->_d_dihedral_coeffs_k4.resize(num_dihedrals_type);

    thrust::copy(fd->_h_dihedral_coeffs_k1,
           fd->_h_dihedral_coeffs_k1 + num_dihedrals_type,
           _device_data->_d_dihedral_coeffs_k1.begin());
    thrust::copy(fd->_h_dihedral_coeffs_k2,
               fd->_h_dihedral_coeffs_k2 + num_dihedrals_type,
               _device_data->_d_dihedral_coeffs_k2.begin());
    thrust::copy(fd->_h_dihedral_coeffs_k3,
               fd->_h_dihedral_coeffs_k3 + num_dihedrals_type,
               _device_data->_d_dihedral_coeffs_k3.begin());
    thrust::copy(fd->_h_dihedral_coeffs_k4,
               fd->_h_dihedral_coeffs_k4 + num_dihedrals_type,
               _device_data->_d_dihedral_coeffs_k4.begin());

  }

  //improper
  std::string improper_type = "null";

  if(*(_structure_info_data->_num_impropers)) {
    improper_type = DataManager::getInstance().getConfigData()->
      Get<std::string>("improper_type", "hyper_parameters", "force_field");
  }
  if (improper_type == "harmonic") {
    _device_data->_d_improper_coeffs_k.resize(num_impropers_type);
    _device_data->_d_improper_coeffs_chi.resize(num_impropers_type);

    thrust::copy(fd->_h_improper_coeffs_k,
               fd->_h_improper_coeffs_k + num_impropers_type,
               _device_data->_d_improper_coeffs_k.begin());
    thrust::copy(fd->_h_improper_coeffs_degree,
                 fd->_h_improper_coeffs_degree + num_impropers_type,
                 _device_data->_d_improper_coeffs_chi.begin());
  }
  else if (improper_type == "cvff") {
    _device_data->_d_improper_coeffs_k.resize(num_impropers_type);
    _device_data->_d_improper_coeffs_d.resize(num_impropers_type);
    _device_data->_d_improper_coeffs_n.resize(num_impropers_type);

    thrust::copy(fd->_h_improper_coeffs_k,
           fd->_h_improper_coeffs_k + num_impropers_type,
           _device_data->_d_improper_coeffs_k.begin());
    thrust::copy(fd->_h_improper_coeffs_d,
               fd->_h_improper_coeffs_d + num_impropers_type,
               _device_data->_d_improper_coeffs_d.begin());
    thrust::copy(fd->_h_improper_coeffs_n,
               fd->_h_improper_coeffs_n + num_impropers_type,
               _device_data->_d_improper_coeffs_n.begin());
  }

  return true;
}

bool CVFFMemoryScheduler::asyncMemoryD2H() { return true; }
