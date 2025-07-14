#pragma once
#include "../../common/device_types.h"
#include "../../common/types.h"
#include "../../data_manager/include/model/box.h"
#include "../force/include/tersoff.h"
namespace op
{
  template <typename DEVICE>
  struct TerSoff
  {
    void operator()(Box box, TersoffParams* params,ShiftFlag shift,
      const rbmd::Real cutmax,
    const rbmd::Id num_atoms,const rbmd::Id nelements,
    const rbmd::Id* atom_id_to_idx,
    const rbmd::Id* atoms_id, const rbmd::Id* atoms_type,const rbmd::Id* map,
    const rbmd::Id* elem3param,
    const rbmd::Id* start_id, const rbmd::Id* end_id,
    const rbmd::Id* id_verletlist, const rbmd::Real* px, const rbmd::Real* py,
    const rbmd::Real* pz,rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,
    rbmd::Real* flat_virial,rbmd::Real* energy);
  };

  template <typename DEVICE>
  struct TerSoffModify
  {
    void operator()(Box box, TersoffParams* params,ShiftFlag shift,
      const rbmd::Real cutmax,
    const rbmd::Id num_atoms,const rbmd::Id nelements,
    const rbmd::Id* atom_id_to_idx,
    const rbmd::Id* atoms_id, const rbmd::Id* atoms_type,const rbmd::Id* map,
    const rbmd::Id* elem3param,
    const rbmd::Id* start_id, const rbmd::Id* end_id,
    const rbmd::Id* id_verletlist, const rbmd::Real* px, const rbmd::Real* py,
    const rbmd::Real* pz,rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,
    rbmd::Real* flat_virial,rbmd::Real* energy,rbmd::Id* d_is_short_neighbor);
  };


// //////////////////////////////////////////////////////
  template <>
  struct TerSoff<device::DEVICE_GPU>
  {
    void operator()( Box box, TersoffParams* params,ShiftFlag shift,
      const rbmd::Real cutmax,
    const rbmd::Id num_atoms, const rbmd::Id nelements,
    const rbmd::Id* atom_id_to_idx,
    const rbmd::Id* atoms_id, const rbmd::Id* atoms_type,const rbmd::Id* map,
    const rbmd::Id* elem3param,
    const rbmd::Id* start_id, const rbmd::Id* end_id,
    const rbmd::Id* id_verletlist, const rbmd::Real* px, const rbmd::Real* py,
    const rbmd::Real* pz,rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,
    rbmd::Real* flat_virial,rbmd::Real* energy);
  };

  template <>
  struct TerSoffModify<device::DEVICE_GPU>
  {
    void operator()(Box box, TersoffParams* params,ShiftFlag shift,
      const rbmd::Real cutmax,
    const rbmd::Id num_atoms,const rbmd::Id nelements,
    const rbmd::Id* atom_id_to_idx,
    const rbmd::Id* atoms_id, const rbmd::Id* atoms_type,const rbmd::Id* map,
    const rbmd::Id* elem3param,
    const rbmd::Id* start_id, const rbmd::Id* end_id,
    const rbmd::Id* id_verletlist, const rbmd::Real* px, const rbmd::Real* py,
    const rbmd::Real* pz,rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,
    rbmd::Real* flat_virial,rbmd::Real* energy,rbmd::Id* d_is_short_neighbor);
  };

}// namespace op