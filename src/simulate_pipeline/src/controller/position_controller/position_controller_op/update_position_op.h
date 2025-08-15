#pragma once
#include "../../data_manager/include/model/box.h"
#include "device_types.h"
#include "types.h"
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
struct UpdatePositionFlagOpvl {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real fmt2v,const rbmd::Real dt,  rbmd::Id test_current_step, Box  box  ,const rbmd::Id* atoms_type,
   const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz,const rbmd::Real* mass,
   rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz,
   rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz,
   rbmd::Real* prev_px,rbmd::Real* prev_py,rbmd::Real* prev_pz,
   rbmd::Id* flag_px,rbmd::Id* flag_py, rbmd::Id* flag_pz);
};

template <typename DEVICE>
struct UpdatePositionFlagOpbm {
  void operator()(const rbmd::Id num_atoms,const rbmd::Real fmt2v, const rbmd::Real dt, rbmd::Id test_current_step,Box  box  , const rbmd::Id* atoms_type,
 const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz,const rbmd::Real* mass,
 rbmd::Real* prev_fx, rbmd::Real* prev_fy, rbmd::Real* prev_fz,
 rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz,
 rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz,
 rbmd::Id* flag_px,rbmd::Id* flag_py, rbmd::Id* flag_pz);
};

template <typename DEVICE>
struct UpdatePositionFlagOp1 {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt, Box box,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* flag_px, rbmd::Id* flag_py,
                  rbmd::Id* flag_pz);
};

template <typename DEVICE>
struct UpdatePositionFlagOp2 {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt, Box box,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* flag_px, rbmd::Id* flag_py,
                  rbmd::Id* flag_pz);
};

template <typename DEVICE>
struct UpdatePositionFlagOp3 {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt, Box box,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* flag_px, rbmd::Id* flag_py,
                  rbmd::Id* flag_pz);
};

template <typename DEVICE>
struct UpdatePositionFlagOp4 {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt, Box box,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* flag_px, rbmd::Id* flag_py,
                  rbmd::Id* flag_pz);
};

// template <typename DEVICE>
// struct UpdatePositionOp {
//   void operator()(const rbmd::Id num_atoms, const rbmd::Real dt, Box box,
//                   const rbmd::Real* vx, const rbmd::Real* vy,
//                   const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
//                   rbmd::Real* pz);
// };

template <>
struct UpdatePositionFlagOp<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt, const rbmd::Real fmt2v, Box box , const rbmd::Id* atoms_type,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, const rbmd::Real* mass,
                  const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz,
                  rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* flag_px, rbmd::Id* flag_py,
                  rbmd::Id* flag_pz);
};

template <>
struct UpdatePositionFlagOp1<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms,const rbmd::Real d1, const rbmd::Real dt, Box box,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* flag_px, rbmd::Id* flag_py,
                  rbmd::Id* flag_pz);
};
template <>
struct UpdatePositionFlagOp2<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms,const rbmd::Real d2, const rbmd::Real dt, Box box,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* flag_px, rbmd::Id* flag_py,
                  rbmd::Id* flag_pz);
};
template <>
struct UpdatePositionFlagOp3<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms,const rbmd::Real d3, const rbmd::Real dt, Box box,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* flag_px, rbmd::Id* flag_py,
                  rbmd::Id* flag_pz);
};
template <>
struct UpdatePositionFlagOp4<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms,const rbmd::Real d4, const rbmd::Real dt, Box box,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* flag_px, rbmd::Id* flag_py,
                  rbmd::Id* flag_pz);
};

template <>
struct UpdatePositionFlagOpvl<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real fmt2v,const rbmd::Real dt, rbmd::Id test_current_step, Box  box  ,const rbmd::Id* atoms_type,
   const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz,const rbmd::Real* mass,
   rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz,
   rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz,
   rbmd::Real* prev_px,rbmd::Real* prev_py,rbmd::Real* prev_pz,
   rbmd::Id* flag_px,rbmd::Id* flag_py, rbmd::Id* flag_pz);
};

template <>
struct UpdatePositionFlagOpbm<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms,const rbmd::Real fmt2v, const rbmd::Real dt, rbmd::Id test_current_step,Box  box  , const rbmd::Id* atoms_type,
 const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz,const rbmd::Real* mass,
 rbmd::Real* prev_fx, rbmd::Real* prev_fy, rbmd::Real* prev_fz,
 rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz,
 rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz,
 rbmd::Id* flag_px,rbmd::Id* flag_py, rbmd::Id* flag_pz);
};

}  // namespace op
