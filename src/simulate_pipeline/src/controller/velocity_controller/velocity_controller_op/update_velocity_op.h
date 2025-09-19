#pragma once
#include "device_types.h"
#include "types.h"
#include "common/device_types.h"


namespace op {
template <typename DEVICE>
struct UpdateVelocityOp {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt,
                  const rbmd::Real fmt2v, const rbmd::Id* atoms_type,
                  const rbmd::Real* mass, const rbmd::Real* fx,
                  const rbmd::Real* fy, const rbmd::Real* fz, rbmd::Real* vx,
                  rbmd::Real* vy, rbmd::Real* vz);
};

template <typename DEVICE>
struct UpdateVelocityOpbm {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real par_a,const rbmd::Real par_b, const rbmd::Real dt, rbmd::Id test_current_step,
                               const rbmd::Real fmt2v,
                               const rbmd::Id* atoms_type,
                               const rbmd::Real* mass, const rbmd::Real* fx,
                               const rbmd::Real* fy, const rbmd::Real* fz,
                               rbmd::Real* prev_fx, rbmd::Real* prev_fy, rbmd::Real* prev_fz,
                               rbmd::Real* pr_prev_fx, rbmd::Real* pr_prev_fy, rbmd::Real* pr_prev_fz,
                               rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz);
};

template <>
struct UpdateVelocityOp<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt,
                  const rbmd::Real fmt2v, const rbmd::Id* atoms_type,
                  const rbmd::Real* mass, const rbmd::Real* fx,
                  const rbmd::Real* fy, const rbmd::Real* fz, rbmd::Real* vx,
                  rbmd::Real* vy, rbmd::Real* vz);
};

template <>
struct UpdateVelocityOpbm<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real par_a,const rbmd::Real par_b, const rbmd::Real dt, rbmd::Id test_current_step,
                               const rbmd::Real fmt2v,
                               const rbmd::Id* atoms_type,
                               const rbmd::Real* mass, const rbmd::Real* fx,
                               const rbmd::Real* fy, const rbmd::Real* fz,
                               rbmd::Real* prev_fx, rbmd::Real* prev_fy, rbmd::Real* prev_fz,
                               rbmd::Real* pr_prev_fx, rbmd::Real* pr_prev_fy, rbmd::Real* pr_prev_fz,
                               rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz);
};

}  // namespace op