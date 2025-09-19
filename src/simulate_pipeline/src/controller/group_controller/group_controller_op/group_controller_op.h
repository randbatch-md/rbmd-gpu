#pragma once
#include "../common/rbmd_define.h"
#include "common/device_types.h"
#include "model/box.h"

namespace op {
template <typename DEVICE>
struct ComputeVCMOp {
  void operator()(const int num_atoms, const int* atoms_type,
      const rbmd::Real* mass,const  rbmd::Real* vx, const  rbmd::Real* vy,
      const rbmd::Real* vz,rbmd::Real* vcm_contrib);
};

template <typename DEVICE>
struct UnwarpPositionOp
{
  void operator()(const rbmd::Id num_atoms, Box  box,
   const rbmd::Real* px, const rbmd::Real* py, const rbmd::Real* pz,
   const rbmd::Id* flag_px, const rbmd::Id* flag_py, const rbmd::Id* flag_pz,
   rbmd::Real* unwarp_px,rbmd::Real* unwarp_py,rbmd::Real* unwarp_pz);
};

template <typename DEVICE>
struct ComputeXCMOp
{
  void operator()(const rbmd::Id num_atoms,
const rbmd::Real* mass,const rbmd::Id* atoms_type,
const rbmd::Real* unwrap_px, const rbmd::Real* unwrap_py,
const rbmd::Real* unwrap_pz,rbmd::Real* xcm_contrib);
};

template <typename DEVICE>
struct ComputeAngMomOp
{
  void operator()( const rbmd::Id num_atoms, Real3 cm,
    const rbmd::Id* d_atoms_type, const rbmd::Real* d_mass,
    const rbmd::Real* d_x, const rbmd::Real* d_y, const rbmd::Real* d_z,
    const rbmd::Real* d_vx, const rbmd::Real* d_vy, const rbmd::Real* d_vz,
    const int* d_image_x, const int* d_image_y, const int* d_image_z,
   Box box,rbmd::Real* d_angmom_contrib);
};

template <typename DEVICE>
struct ComputeInertiaOp
{
  void operator()( const rbmd::Id num_atoms, Real3 cm,
    const rbmd::Id* d_atoms_type, const rbmd::Real* d_mass,
    const rbmd::Real* d_x, const rbmd::Real* d_y, const rbmd::Real* d_z,
    const int* d_image_x, const int* d_image_y, const int* d_image_z,
   Box box,rbmd::Real* d_inertia_contrib);
};


//////////////////////////////////////////
template <>
struct ComputeVCMOp<device::DEVICE_GPU> {
  void operator()(const int num_atoms, const int* atoms_type,
      const rbmd::Real* mass,const  rbmd::Real* vx, const  rbmd::Real* vy,
      const rbmd::Real* vz,rbmd::Real* vcm_contrib);
};

template <>
struct UnwarpPositionOp<device::DEVICE_GPU>
{
  void operator()(const rbmd::Id num_atoms, Box box, const rbmd::Real* px,
                  const rbmd::Real* py, const rbmd::Real* pz,
                  const rbmd::Id* flag_px, const rbmd::Id* flag_py,
                  const rbmd::Id* flag_pz, rbmd::Real* unwarp_px,
                  rbmd::Real* unwarp_py, rbmd::Real* unwarp_pz);
};

template <>
struct ComputeXCMOp<device::DEVICE_GPU>
{
  void operator()(const rbmd::Id num_atoms,
const rbmd::Real* mass,const rbmd::Id* atoms_type,
const rbmd::Real* unwrap_px, const rbmd::Real* unwrap_py,
const rbmd::Real* unwrap_pz,rbmd::Real* xcm_contrib);
};

template <>
struct ComputeAngMomOp<device::DEVICE_GPU>
{
  void operator()( const rbmd::Id num_atoms, Real3 cm,
    const rbmd::Id* d_atoms_type, const rbmd::Real* d_mass,
    const rbmd::Real* d_x, const rbmd::Real* d_y, const rbmd::Real* d_z,
    const rbmd::Real* d_vx, const rbmd::Real* d_vy, const rbmd::Real* d_vz,
    const int* d_image_x, const int* d_image_y, const int* d_image_z,
    Box box,rbmd::Real* d_angmom_contrib);
};

template <>
struct ComputeInertiaOp<device::DEVICE_GPU>
{
  void operator()( const rbmd::Id num_atoms, Real3 cm,
    const rbmd::Id* d_atoms_type, const rbmd::Real* d_mass,
    const rbmd::Real* d_x, const rbmd::Real* d_y, const rbmd::Real* d_z,
    const int* d_image_x, const int* d_image_y, const int* d_image_z,
   Box box,rbmd::Real* d_inertia_contrib);
};
}