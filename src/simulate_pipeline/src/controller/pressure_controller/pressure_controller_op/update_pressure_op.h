#pragma once
#include "device_types.h"
#include "types.h"

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

}