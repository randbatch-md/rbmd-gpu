#pragma once
#include "../common/rbmd_define.h"
#include "../common/types.h"
#include "box.h"
#include "../common/erf_table.h"

class DeviceData {
 public:
  ///< struct data on device>
  /// position
  thrust::device_vector<rbmd::Real> _d_px;
  thrust::device_vector<rbmd::Real> _d_py;
  thrust::device_vector<rbmd::Real> _d_pz;

  thrust::device_vector<rbmd::Real> _d_shake_px;
  thrust::device_vector<rbmd::Real> _d_shake_py;
  thrust::device_vector<rbmd::Real> _d_shake_pz;

  /// atoms id
  thrust::device_vector<rbmd::Id> _d_atoms_id;

  /// atoms type
  thrust::device_vector<rbmd::Id> _d_atoms_type;

  /// molecular id
  thrust::device_vector<rbmd::Id> _d_molecular_id;

  /// atoms flag
  thrust::device_vector<rbmd::Id> _d_flagX;
  thrust::device_vector<rbmd::Id> _d_flagY;
  thrust::device_vector<rbmd::Id> _d_flagZ;

  /// velocity
  thrust::device_vector<rbmd::Real> _d_vx;
  thrust::device_vector<rbmd::Real> _d_vy;
  thrust::device_vector<rbmd::Real> _d_vz;
  thrust::device_vector<rbmd::Real> _d_shake_vx;
  thrust::device_vector<rbmd::Real> _d_shake_vy;
  thrust::device_vector<rbmd::Real> _d_shake_vz;

  /// bond
  thrust::device_vector<rbmd::Id> _d_bond_type;
  thrust::device_vector<rbmd::Id> _d_bond_id0;
  thrust::device_vector<rbmd::Id> _d_bond_id1;

  thrust::device_vector<rbmd::Real> _d_special_weights;
  thrust::device_vector<rbmd::Id> _d_special_ids;
  thrust::device_vector<rbmd::Id> _d_special_count;
  thrust::device_vector<rbmd::Id> _d_special_offsets;

  thrust::device_vector<rbmd::Id> _d_atoms_vec;
  thrust::device_vector<rbmd::Id> _d_atoms_count;
  thrust::device_vector<rbmd::Id> _d_atoms_offset;

  /// angle
  thrust::device_vector<rbmd::Id> _d_angle_type;
  thrust::device_vector<rbmd::Id> _d_angle_id0;
  thrust::device_vector<rbmd::Id> _d_angle_id1;
  thrust::device_vector<rbmd::Id> _d_angle_id2;
  thrust::device_vector<Id3> _d_angle_id_vec;

  /// dihedral
  thrust::device_vector<rbmd::Id> _d_dihedral_type;
  thrust::device_vector<rbmd::Id> _d_dihedral_id0;
  thrust::device_vector<rbmd::Id> _d_dihedral_id1;
  thrust::device_vector<rbmd::Id> _d_dihedral_id2;
  thrust::device_vector<rbmd::Id> _d_dihedral_id3;

  /// Processing data
  // thrust::device_vector<rbmd::Id> _d_special_source_array;
  // thrust::device_vector<rbmd::Id> _d_special_offsets_array;


  ///< force on device>
  /// mass
  thrust::device_vector<rbmd::Real> _d_mass;

  /// eps
  thrust::device_vector<rbmd::Real> _d_eps;

  /// sigma
  thrust::device_vector<rbmd::Real> _d_sigma;

  /// bond
  thrust::device_vector<rbmd::Real> _d_bond_coeffs_k;
  thrust::device_vector<rbmd::Real> _d_bond_coeffs_equilibrium;

  /// angle
  thrust::device_vector<rbmd::Real> _d_angle_coeffs_k;
  thrust::device_vector<rbmd::Real> _d_angle_coeffs_equilibrium;

  /// dihedral
  thrust::device_vector<rbmd::Real> _d_dihedral_coeffs_k;
  thrust::device_vector<rbmd::Id> _d_dihedral_coeffs_sign;
  thrust::device_vector<rbmd::Id> _d_dihedral_coeffs_multiplicity;

  /// F(ρ) on host
  thrust::device_vector<rbmd::Real> _d_frho;

  /// ρ(r) on host
  thrust::device_vector<rbmd::Real> _d_rhor;

  /// ϕ(r) on host
  thrust::device_vector<rbmd::Real> _d_z2r;

  ///< force on device>
  thrust::device_vector<rbmd::Real> _d_fx;
  thrust::device_vector<rbmd::Real> _d_fy;
  thrust::device_vector<rbmd::Real> _d_fz;
  //lj_cut + coul_cut
  thrust::device_vector<rbmd::Real> _d_force_ljcoul_x;
  thrust::device_vector<rbmd::Real> _d_force_ljcoul_y;
  thrust::device_vector<rbmd::Real> _d_force_ljcoul_z;
  //kspace
  thrust::device_vector<rbmd::Real> _d_force_kspace_x;
  thrust::device_vector<rbmd::Real> _d_force_kspace_y;
  thrust::device_vector<rbmd::Real> _d_force_kspace_z;
  //special_coul
  thrust::device_vector<rbmd::Real> _d_force_specialcoul_x;
  thrust::device_vector<rbmd::Real> _d_force_specialcoul_y;
  thrust::device_vector<rbmd::Real> _d_force_specialcoul_z;
  //bond
  thrust::device_vector<rbmd::Real> _d_force_bond_x;
  thrust::device_vector<rbmd::Real> _d_force_bond_y;
  thrust::device_vector<rbmd::Real> _d_force_bond_z;

  thrust::device_vector<rbmd::Id> _d_temp_atom_ids;
  thrust::device_vector<rbmd::Real> _d_temp_forces_bondx;
  thrust::device_vector<rbmd::Real> _d_temp_forces_bondy;
  thrust::device_vector<rbmd::Real> _d_temp_forces_bondz;
  //angle
  thrust::device_vector<rbmd::Real> _d_force_angle_x;
  thrust::device_vector<rbmd::Real> _d_force_angle_y;
  thrust::device_vector<rbmd::Real> _d_force_angle_z;
  //dihedral
  thrust::device_vector<rbmd::Real> _d_force_dihedral_x;
  thrust::device_vector<rbmd::Real> _d_force_dihedral_y;
  thrust::device_vector<rbmd::Real> _d_force_dihedral_z;

  //virial
  thrust::device_vector<rbmd::Real> _d_flat_virial;
  thrust::device_vector<rbmd::Real> _d_flat_virial6;
  thrust::device_vector<rbmd::Real> _d_flat_virial_lj;
  thrust::device_vector<rbmd::Real> _d_flat_virial_specialcoul;
  thrust::device_vector<rbmd::Real> _d_flat_virial_kspace;

  thrust::device_vector<rbmd::Real> _d_flat_virial_bond;
  thrust::device_vector<rbmd::Real> _d_flat_virial_angle;
  thrust::device_vector<rbmd::Real> _d_flat_virial_dihedral;

  thrust::device_vector<rbmd::Real> _d_flat_virial_bond_atom;
  thrust::device_vector<rbmd::Real> _d_flat_virial_angle_atom;

  thrust::device_vector<rbmd::Real> _d_virial;
  thrust::device_vector<rbmd::Real> _d_virial_lj;
  thrust::device_vector<rbmd::Real> _d_virial_specialcoul;
  thrust::device_vector<rbmd::Real> _d_virial_kspace;
  thrust::device_vector<rbmd::Real> _d_virial_bond;
  thrust::device_vector<rbmd::Real> _d_virial_angle;
  thrust::device_vector<rbmd::Real> _d_virial_dihedral;
  thrust::device_vector<rbmd::Real> _d_energy_dihedral;

  //charge
  thrust::device_vector<rbmd::Real> _d_charge;

  // ERFTable in device
  ERFTable* _d_erf_table;

  void unload()
  {
	  CHECK_RUNTIME(FREE(_d_erf_table));
  }
};
