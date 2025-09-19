#pragma once
#include "common/device_types.h"
#include  "../common/rbmd_define.h"
#include "../data_manager/include/model/box.h"
namespace op {

template <typename DEVICE>
struct ComputeKineticEnergyOp
{
  void operator()(const rbmd::Id num_atoms,rbmd::Id* atoms_type,
    const rbmd::Real* mass,const rbmd::Real* vx, const rbmd::Real* vy,
    const rbmd::Real* vz,rbmd::Real* d_ke_contrib);
};

template <typename DEVICE>
struct ZeroLinearMomentumOp
{
  void operator()(const rbmd::Id num_atoms,Real3 vcm,
  bool x_flag, bool y_flag, bool z_flag,rbmd::Real* vx,
  rbmd::Real* vy,rbmd::Real* vz);
};



template <typename DEVICE>
struct ZeroAngularMomentumOp
{
  void operator()( const rbmd::Id num_atoms, const rbmd::Real* xcm, const rbmd::Real* omega,
    const rbmd::Real* d_x, const rbmd::Real* d_y, const rbmd::Real* d_z,
    const int* d_image_x, const int* d_image_y, const int* d_image_z,
    const rbmd::Real* d_box_prd,
    rbmd::Real* d_vx, rbmd::Real* d_vy, rbmd::Real* d_vz);
};


/////////////////////////
template <>
struct ComputeKineticEnergyOp<device::DEVICE_GPU>
{
  void operator()(const rbmd::Id num_atoms,rbmd::Id* atoms_type,
    const rbmd::Real* mass,const rbmd::Real* vx, const rbmd::Real* vy,
    const rbmd::Real* vz,rbmd::Real* d_ke_contrib);
};

template <>
struct ZeroLinearMomentumOp<device::DEVICE_GPU>
{
  void operator()(const rbmd::Id num_atoms,Real3 vcm,
  bool x_flag, bool y_flag, bool z_flag,rbmd::Real* vx,
  rbmd::Real* vy,rbmd::Real* vz);
};



template <>
struct ZeroAngularMomentumOp<device::DEVICE_GPU>
{
  void operator()( const rbmd::Id num_atoms, Real3 xcm, const rbmd::Real* omega,
    const rbmd::Real* px, const rbmd::Real* py, const rbmd::Real* pz,
    const int* d_image_x, const int* d_image_y, const int* d_image_z,
    Box box, rbmd::Real* d_vx, rbmd::Real* d_vy, rbmd::Real* d_vz);
};

}