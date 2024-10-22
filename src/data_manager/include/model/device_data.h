#pragma once
#include <thrust/device_vector.h>

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

  thrust::device_vector<rbmd::Real> _d_force_ljcoul_x;
  thrust::device_vector<rbmd::Real> _d_force_ljcoul_y;
  thrust::device_vector<rbmd::Real> _d_force_ljcoul_z;

  thrust::device_vector<rbmd::Real> _d_force_ewald_x;
  thrust::device_vector<rbmd::Real> _d_force_ewald_y;
  thrust::device_vector<rbmd::Real> _d_force_ewald_z;

  thrust::device_vector<rbmd::Real> _d_force_specialcoul_x;
  thrust::device_vector<rbmd::Real> _d_force_specialcoul_y;
  thrust::device_vector<rbmd::Real> _d_force_specialcoul_z;
  //
  thrust::device_vector<rbmd::Real> _d_force_bond_x;
  thrust::device_vector<rbmd::Real> _d_force_bond_y;
  thrust::device_vector<rbmd::Real> _d_force_bond_z;

  thrust::device_vector<rbmd::Real> _d_force_angle_x;
  thrust::device_vector<rbmd::Real> _d_force_angle_y;
  thrust::device_vector<rbmd::Real> _d_force_angle_z;



  thrust::device_vector<rbmd::Real> _d_virial_xx;
  thrust::device_vector<rbmd::Real> _d_virial_yy;
  thrust::device_vector<rbmd::Real> _d_virial_zz;
  thrust::device_vector<rbmd::Real> _d_virial_xy;
  thrust::device_vector<rbmd::Real> _d_virial_xz;
  thrust::device_vector<rbmd::Real> _d_virial_yz;
  //
  thrust::device_vector<rbmd::Real> _d_charge;

  // box in device
  Box* _d_box;
  ERFTable* _d_erf_table;
  ~DeviceData() { CHECK_RUNTIME(FREE(_d_box));CHECK_RUNTIME(FREE(_d_erf_table)); }
};
