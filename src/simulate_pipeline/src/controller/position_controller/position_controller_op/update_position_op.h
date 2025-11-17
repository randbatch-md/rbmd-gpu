#pragma once
#include "../../data_manager/include/model/box.h"
//#include "types.h"
#include "../common/types.h"
#include "common/device_types.h"

namespace op {
template <typename DEVICE>
struct UpdatePositionFlagOp {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt, Box box,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* flag_px, rbmd::Id* flag_py,
                  rbmd::Id* flag_pz);
};

template <typename DEVICE>
struct PBCOp {
  void operator()(const rbmd::Id num_atoms,  Box box,
    rbmd::Real* px, rbmd::Real* py,rbmd::Real* pz, rbmd::Id* flag_px,
    rbmd::Id* flag_py,rbmd::Id* flag_pz);
};

template <typename DEVICE>
struct UpdatePositionOp {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt, Box box,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz);
};

template <typename DEVICE>
struct UpdatePositionFlagOpbm {
  void operator()(const rbmd::Id num_atoms,const rbmd::Real fmt2v,
     const rbmd::Real par_a,const rbmd::Real par_b, const rbmd::Real dt,
     rbmd::Id test_current_step,Box  box  , const rbmd::Id* atoms_type,
     const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz,
     const rbmd::Real* mass,rbmd::Real* prev_fx, rbmd::Real* prev_fy, rbmd::Real* prev_fz,
     rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz,
     rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz,
     rbmd::Id* flag_px,rbmd::Id* flag_py, rbmd::Id* flag_pz);
};

template <>
struct UpdatePositionFlagOp<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt, Box box,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* flag_px, rbmd::Id* flag_py,
                  rbmd::Id* flag_pz);
};

template <>
struct PBCOp<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms,  Box box,
    rbmd::Real* px, rbmd::Real* py,rbmd::Real* pz, rbmd::Id* flag_px,
    rbmd::Id* flag_py,rbmd::Id* flag_pz);
};

template <>
struct UpdatePositionOp<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt, Box box,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz);
};

template <>
struct UpdatePositionFlagOpbm<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms,const rbmd::Real fmt2v,
     const rbmd::Real par_a,const rbmd::Real par_b, const rbmd::Real dt,
     rbmd::Id test_current_step,Box  box  , const rbmd::Id* atoms_type,
     const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz,const rbmd::Real* mass,
     rbmd::Real* prev_fx, rbmd::Real* prev_fy, rbmd::Real* prev_fz,
     rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz,
     rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz,
     rbmd::Id* flag_px,rbmd::Id* flag_py, rbmd::Id* flag_pz);
};

}  // namespace op
