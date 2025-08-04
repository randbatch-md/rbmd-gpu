#pragma once
#include "../../data_manager/include/model/box.h"
#include "device_types.h"
#include "types.h"
#include "common/device_types.h"

namespace op {
//leapfrog
template <typename DEVICE>
struct UpdatePositionFlagOp {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt, Box box,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* flag_px, rbmd::Id* flag_py,
                  rbmd::Id* flag_pz);
};

//vv
template <typename DEVICE>
struct UpdatePositionFlagOpvv {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt, Box box,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* flag_px, rbmd::Id* flag_py,
                  rbmd::Id* flag_pz);
};

//PRK
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

//未调用
template <typename DEVICE>
struct UpdatePositionOp {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt, Box box,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz);
};

//leapfrog
template <>
struct UpdatePositionFlagOp<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt,
                  const rbmd::Real fmt2v, const rbmd::Id* atoms_type, const rbmd::Real* mass,
                  Box box,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* flag_px, rbmd::Id* flag_py,
                  rbmd::Id* flag_pz,
                  const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz);
};

//vv
template <>
struct UpdatePositionFlagOpvv<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt,
                  const rbmd::Real fmt2v, const rbmd::Id* atoms_type, const rbmd::Real* mass,
                  Box box,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* flag_px, rbmd::Id* flag_py,
                  rbmd::Id* flag_pz,
                  const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz);
};

// PRK
template <>
struct UpdatePositionFlagOp1<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real d1, const rbmd::Real dt,
                  const rbmd::Real fmt2v, const rbmd::Id* atoms_type, const rbmd::Real* mass,
                  Box box,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* flag_px, rbmd::Id* flag_py,
                  rbmd::Id* flag_pz,
                  const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz);
};

template <>
struct UpdatePositionFlagOp2<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real d2, const rbmd::Real dt,
                  const rbmd::Real fmt2v, const rbmd::Id* atoms_type, const rbmd::Real* mass,
                  Box box,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* flag_px, rbmd::Id* flag_py,
                  rbmd::Id* flag_pz,
                  const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz);
};

template <>
struct UpdatePositionFlagOp3<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real d3, const rbmd::Real dt,
                  const rbmd::Real fmt2v, const rbmd::Id* atoms_type, const rbmd::Real* mass,
                  Box box,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* flag_px, rbmd::Id* flag_py,
                  rbmd::Id* flag_pz,
                  const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz);
};

template <>
struct UpdatePositionFlagOp4<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real d4, const rbmd::Real dt,
                  const rbmd::Real fmt2v, const rbmd::Id* atoms_type, const rbmd::Real* mass,
                  Box box,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* flag_px, rbmd::Id* flag_py,
                  rbmd::Id* flag_pz,
                  const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz);
};

//未调用
template <>
struct UpdatePositionOp<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt, Box box,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz);
};

}  // namespace op
