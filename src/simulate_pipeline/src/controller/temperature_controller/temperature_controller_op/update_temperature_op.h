#pragma once
#include "device_types.h"
#include "types.h"

namespace op {
template <typename DEVICE>
struct ComputeTemperatureOp {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real mvv2e,
                  const rbmd::Id* atoms_type, const rbmd::Real* mass,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, rbmd::Real* temp_contrib);
};

template <typename DEVICE>
struct UpdataVelocityRescaleOp {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real coeff_rescale,
                  rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz);
};

template <typename DEVICE>
struct UpdataVelocityNoseHooverOp {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt,
                  const rbmd::Real fmt2v, const rbmd::Real nosehooverxi,
                  const rbmd::Id* atoms_type, const rbmd::Real* mass,
                  const rbmd::Real* fx, const rbmd::Real* fy,
                  const rbmd::Real* fz, rbmd::Real* vx, rbmd::Real* vy,
                  rbmd::Real* vz);
};

template <typename DEVICE>
struct UpdataVelocityBerendsenOp {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real coeff_Berendsen,
                  rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz);
};

template <>
struct ComputeTemperatureOp<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real mvv2e,
                  const rbmd::Id* atoms_type, const rbmd::Real* mass,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, rbmd::Real* temp_contrib);
};

template <>
struct UpdataVelocityRescaleOp<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real coeff_rescale,
                  rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz);
};

template <>
struct UpdataVelocityNoseHooverOp<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt,
                  const rbmd::Real fmt2v, const rbmd::Real nosehooverxi,
                  const rbmd::Id* atoms_type, const rbmd::Real* mass,
                  const rbmd::Real* fx, const rbmd::Real* fy,
                  const rbmd::Real* fz, rbmd::Real* vx, rbmd::Real* vy,
                  rbmd::Real* vz);
};

template <>
struct UpdataVelocityBerendsenOp<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real coeff_Berendsen,
                  rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz);
};

}  // namespace op