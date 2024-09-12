#pragma once
#include "device_types.h"
#include "types.h"
namespace op {
template <typename DEVICE>
struct UpdateVelocityOp {
  void operator()(const rbmd::Id num_atoms, 
                  const rbmd::Real dt,
                  const rbmd::Real fmt2v, 
                  const rbmd::Real* mass,
                  const rbmd::Real* fx, 
                  const rbmd::Real* fy,
                  const rbmd::Real* fz, 
                  rbmd::Real* vx, 
                  rbmd::Real* vy,
                  rbmd::Real* vz);
};

template <>
struct UpdateVelocityOp<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms, 
                  const rbmd::Real dt,
                  const rbmd::Real fmt2v, 
                  const rbmd::Real* mass,
                  const rbmd::Real* fx, 
                  const rbmd::Real* fy,
                  const rbmd::Real* fz, 
                  rbmd::Real* vx, 
                  rbmd::Real* vy,
                  rbmd::Real* vz);
};
}  // namespace op