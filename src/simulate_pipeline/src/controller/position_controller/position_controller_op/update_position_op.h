#pragma once
#include "../../data_manager/include/model/box.h"
#include "device_types.h"
#include "types.h"

namespace op {
template <typename DEVICE>
struct UpdatePositionFlagOp {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt, Box* box,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* flag_px, rbmd::Id* flag_py,
                  rbmd::Id* flag_pz);
};

template <typename DEVICE>
struct UpdatePositionOp {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt, Box* box,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz);
};

template <>
struct UpdatePositionFlagOp<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt, Box* box,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* flag_px, rbmd::Id* flag_py,
                  rbmd::Id* flag_pz);
};

template <>
struct UpdatePositionOp<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms, const rbmd::Real dt, Box* box,
                  const rbmd::Real* vx, const rbmd::Real* vy,
                  const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz);
};
}  // namespace op
