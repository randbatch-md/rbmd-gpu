#pragma once

#include <thrust/host_vector.h>

#include <fstream>

#include "common/types.h"
#include "force_field_data.h"

class EAMForceFieldData : public ForceFieldData {
 public:
  bool checkForceField() const override { return true; }


  rbmd::Real _cut_off;

  /// F(ρ) on host
  rbmd::Id _nrho;
  rbmd::Real _drho;
  //rbmd::Real* _h_frho;
  thrust::host_vector<rbmd::Real> _h_frho;


  /// ρ(r) on host
  // rbmd::Real* _h_rhor;
  thrust::host_vector<rbmd::Real> _h_rhor;

  /// ϕ(r) on host
  rbmd::Id _nr;
  rbmd::Real _dr;
  // rbmd::Real* _h_zr;
  thrust::host_vector<rbmd::Real> _h_zr;

};
