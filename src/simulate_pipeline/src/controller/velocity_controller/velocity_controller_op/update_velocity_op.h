#pragma once
#include "device_types.h"
#include "types.h"
#include "common/device_types.h"


namespace op {
// leapfrog
template <typename DEVICE>
struct UpdateVelocityOp {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt,
                  const rbmd::Real fmt2v, const rbmd::Id* atoms_type,
                  const rbmd::Real* mass, const rbmd::Real* fx,
                  const rbmd::Real* fy, const rbmd::Real* fz, rbmd::Real* vx,
                  rbmd::Real* vy, rbmd::Real* vz);
};

// vv
template <typename DEVICE>
struct UpdateVelocityOpvv {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt,
                  const rbmd::Real fmt2v, const rbmd::Id* atoms_type,
                  const rbmd::Real* mass, const rbmd::Real* fx,
                  const rbmd::Real* fy, const rbmd::Real* fz, rbmd::Real* vx,
                  rbmd::Real* vy, rbmd::Real* vz);
};

// PRK
template <typename DEVICE>
struct UpdateVelocityOp1 {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real c1, const rbmd::Real dt,
                  const rbmd::Real fmt2v, const rbmd::Id* atoms_type,
                  const rbmd::Real* mass, const rbmd::Real* fx,
                  const rbmd::Real* fy, const rbmd::Real* fz, rbmd::Real* vx,
                  rbmd::Real* vy, rbmd::Real* vz);
};
template <typename DEVICE>
struct UpdateVelocityOp2 {
  void operator()(const rbmd::Id num_atoms,const rbmd::Real c2, const rbmd::Real dt,
                  const rbmd::Real fmt2v, const rbmd::Id* atoms_type,
                  const rbmd::Real* mass, const rbmd::Real* fx,
                  const rbmd::Real* fy, const rbmd::Real* fz, rbmd::Real* vx,
                  rbmd::Real* vy, rbmd::Real* vz);
};
template <typename DEVICE>
struct UpdateVelocityOp3 {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt,
                  const rbmd::Real fmt2v, const rbmd::Id* atoms_type,
                  const rbmd::Real* mass, const rbmd::Real* fx,
                  const rbmd::Real* fy, const rbmd::Real* fz, rbmd::Real* vx,
                  rbmd::Real* vy, rbmd::Real* vz);
};
template <typename DEVICE>
struct UpdateVelocityOp4 {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt,
                  const rbmd::Real fmt2v, const rbmd::Id* atoms_type,
                  const rbmd::Real* mass, const rbmd::Real* fx,
                  const rbmd::Real* fy, const rbmd::Real* fz, rbmd::Real* vx,
                  rbmd::Real* vy, rbmd::Real* vz);
};

//beeman
template <typename DEVICE>
struct UpdateVelocityOpBeeman {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt,
                  const rbmd::Real fmt2v, const rbmd::Id* atoms_type,
                  const rbmd::Real* mass, const rbmd::Real* fx,
                  const rbmd::Real* fy, const rbmd::Real* fz, rbmd::Real* vx,
                  rbmd::Real* vy, rbmd::Real* vz);
};

/// GPU
//leapfrog
template <>
struct UpdateVelocityOp<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt,
                  const rbmd::Real fmt2v, const rbmd::Id* atoms_type,
                  const rbmd::Real* mass, const rbmd::Real* fx,
                  const rbmd::Real* fy, const rbmd::Real* fz, rbmd::Real* vx,
                  rbmd::Real* vy, rbmd::Real* vz);
};

//vv
template <>
struct UpdateVelocityOpvv<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt,
                  const rbmd::Real fmt2v, const rbmd::Id* atoms_type, const rbmd::Real* mass,
                  const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz,
                  const rbmd::Real* fx_prev, const rbmd::Real* fy_prev, const rbmd::Real* fz_prev,
                  rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz);
};

//PRK
template <>
struct UpdateVelocityOp1<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real c1, const rbmd::Real dt,
                  const rbmd::Real fmt2v, const rbmd::Id* atoms_type,
                  const rbmd::Real* mass, const rbmd::Real* fx,
                  const rbmd::Real* fy, const rbmd::Real* fz, rbmd::Real* vx,
                  rbmd::Real* vy, rbmd::Real* vz);
};

template <>
struct UpdateVelocityOp2<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real c2, const rbmd::Real dt,
                  const rbmd::Real fmt2v, const rbmd::Id* atoms_type,
                  const rbmd::Real* mass, const rbmd::Real* fx,
                  const rbmd::Real* fy, const rbmd::Real* fz, rbmd::Real* vx,
                  rbmd::Real* vy, rbmd::Real* vz);
};

template <>
struct UpdateVelocityOp3<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real c3, const rbmd::Real dt,
                  const rbmd::Real fmt2v, const rbmd::Id* atoms_type,
                  const rbmd::Real* mass, const rbmd::Real* fx,
                  const rbmd::Real* fy, const rbmd::Real* fz, rbmd::Real* vx,
                  rbmd::Real* vy, rbmd::Real* vz);
};

template <>
struct UpdateVelocityOp4<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real c4, const rbmd::Real dt,
                  const rbmd::Real fmt2v, const rbmd::Id* atoms_type,
                  const rbmd::Real* mass, const rbmd::Real* fx,
                  const rbmd::Real* fy, const rbmd::Real* fz, rbmd::Real* vx,
                  rbmd::Real* vy, rbmd::Real* vz);
};

//beeman
template <>
struct UpdateVelocityOpBeeman<device::DEVICE_GPU> {
  void operator()(
      const rbmd::Id num_atoms, const rbmd::Real dt, rbmd::Id test_current_step,const rbmd::Real fmt2v,
      const rbmd::Id* atoms_type, const rbmd::Real* mass,
      const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz,
      rbmd::Real* f_pre1_x, rbmd::Real* f_pre1_y, rbmd::Real* f_pre1_z,
      rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz);
};

}  // namespace op