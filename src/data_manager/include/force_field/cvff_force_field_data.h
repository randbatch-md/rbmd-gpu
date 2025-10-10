#pragma once

#include "common/types.h"
#include "force_field_data.h"

class CVFFForceFieldData : public ForceFieldData {
 public:
  bool checkForceField() const override { return true; }

  /// mass
  //rbmd::Real* _h_mass;

  /// eps
  rbmd::Real* _h_eps;

  /// sigma
  rbmd::Real* _h_sigma;

  /// bond
  rbmd::Real* _h_bond_coeffs_k;
  rbmd::Real* _h_bond_coeffs_equilibrium;

  /// angle
  rbmd::Real* _h_angle_coeffs_k;
  rbmd::Real* _h_angle_coeffs_equilibrium;

  /// dihedral
  rbmd::Real* _h_dihedral_coeffs_k;
  rbmd::Id* _h_dihedral_coeffs_sign;
  rbmd::Id* _h_dihedral_coeffs_multiplicity;

  rbmd::Id*  _h_nterms;
  rbmd::Id* _h_fourier_offsets;
  rbmd::Real* _h_fourier_cos_shift;
  rbmd::Real* _h_fourier_sin_shift;

  rbmd::Real* _h_dihedral_coeffs_k1;
  rbmd::Real* _h_dihedral_coeffs_k2;
  rbmd::Real* _h_dihedral_coeffs_k3;
  rbmd::Real* _h_dihedral_coeffs_k4;

  //improper
  rbmd::Real* _h_improper_coeffs_k;
  rbmd::Real* _h_improper_coeffs_degree;

  rbmd::Id* _h_improper_coeffs_d;
  rbmd::Id* _h_improper_coeffs_n;
};
