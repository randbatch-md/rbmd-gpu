#pragma once

#include "types.h"
#include  "../common/device_types.h"
#include "model/box.h"

namespace op {

template <typename DEVICE>
struct X2LamdaOp
{
  void operator()(Box box,
                  const rbmd::Id num_atoms,
                   rbmd::Real* px,
                   rbmd::Real* py,
                   rbmd::Real* pz);
};

template <typename DEVICE>
struct Lamda2XOp
{
  void operator()(Box box,
                  const rbmd::Id num_atoms,
                   rbmd::Real* px,
                   rbmd::Real* py,
                   rbmd::Real* pz);
};

template <typename DEVICE>
struct UpdataVelocityRescalePressureOp
{
  void operator()(const rbmd::Id num_atoms, const Real3 factor,
                  rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz);
};

//////////////

template <>
struct X2LamdaOp<device::DEVICE_GPU>
{
  void operator()(Box box,
                  const rbmd::Id num_atoms,
                   rbmd::Real* px,
                   rbmd::Real* py,
                   rbmd::Real* pz);
};

template <>
struct Lamda2XOp<device::DEVICE_GPU>
{
  void operator()(Box box,
                  const rbmd::Id num_atoms,
                   rbmd::Real* px,
                   rbmd::Real* py,
                   rbmd::Real* pz);
};

template <>
struct UpdataVelocityRescalePressureOp<device::DEVICE_GPU> {
  void operator()(const rbmd::Id num_atoms, const Real3 factor,
                  rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz);
};

}