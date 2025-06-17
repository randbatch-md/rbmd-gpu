#pragma once
#include "../../common/device_types.h"
#include "../../common/types.h"
#include "../../data_manager/include/model/box.h"
#include "../../force/include/eam.h"

namespace op {

//Verlet:: Fp Verlet
template <typename DEVICE>
struct ComputeFpVerlet
{
  void operator()( Box box,
        EAMParameters eam_paras,
        const rbmd::Real cut_off,
        const rbmd::Id num_atoms,
        const rbmd::Id* atoms_type,
        const rbmd::Id* atoms_id,
        const rbmd::Id* start_id,
        const rbmd::Id* end_id,
        const rbmd::Id* id_verletlist,
        const Real7*  rhor_spline,
        const Real7* frho_spline,
        const rbmd::Real* px,
        const rbmd::Real* py,
        const rbmd::Real* pz,
        rbmd::Real* eam_fp,
        rbmd::Real* energy_embedding);
};

//Verlet:: EAMForceVerletOp
template <typename DEVICE>
struct ComputeEAMForceVerlet
{
  void operator()( Box box,
        EAMParameters eam_paras,
        const rbmd::Real cut_off,
        const rbmd::Id num_atoms,
        const rbmd::Id* atoms_type,
        const rbmd::Id* atoms_id,
        const rbmd::Id* start_id,
        const rbmd::Id* end_id,
        const rbmd::Id* id_verletlist,
        const Real7*  rhor_spline,
        const Real7* frho_spline,
        const Real7*  z2r_spline,
        const rbmd::Real* px,
        const rbmd::Real* py,
        const rbmd::Real* pz,
        rbmd::Real* eam_rho,
        rbmd::Real* eam_fp,
        rbmd::Real* fx,
        rbmd::Real* fy,
        rbmd::Real* fz,
        rbmd::Real* energy_embedding,
        rbmd::Real* energy_pair);
};

//RBL :: EAMRhoRBL
template <typename DEVICE>
struct ComputeRhoRBL
{
  void operator()( Box box,
        EAMParameters eam_paras,
        const rbmd::Real rs,
        const rbmd::Real rc,
        const rbmd::Id num_atoms,
        const rbmd::Id neighbor_sample_num,
        const rbmd::Id pice_num,
        const rbmd::Id* atoms_type,
        const rbmd::Id* atoms_id,
        const rbmd::Id* start_id,
        const rbmd::Id* end_id,
        const rbmd::Id* id_verletlist,
        const rbmd::Id* id_random_neighbor,
        const rbmd::Id* random_neighbor_num,
        const Real7*  rhor_spline,
        const rbmd::Real* px,
        const rbmd::Real* py,
        const rbmd::Real* pz,
        rbmd::Real* eam_rho);
};


template <typename DEVICE>
struct ComputeFpRBL
{
  void operator()( Box box,
        EAMParameters eam_paras,
        const rbmd::Real rs,
        const rbmd::Real rc,
        const rbmd::Id num_atoms,
        const rbmd::Id neighbor_sample_num,
        const rbmd::Id pice_num,
        const rbmd::Id* atoms_type,
        const rbmd::Id* atoms_id,
        const rbmd::Id* start_id,
        const rbmd::Id* end_id,
        const rbmd::Id* id_verletlist,
        const rbmd::Id* id_random_neighbor,
        const rbmd::Id* random_neighbor_num,
        const Real7*  rhor_spline,
        const Real7* frho_spline,
        const rbmd::Real* px,
        const rbmd::Real* py,
        const rbmd::Real* pz,
        rbmd::Real* eam_fp);
};

template <typename DEVICE>
struct ComputeEAMForceRBL
{
  void operator()( Box box,
        EAMParameters eam_paras,
        const rbmd::Real rs,
        const rbmd::Real rc,
        const rbmd::Id num_atoms,
        const rbmd::Id neighbor_sample_num,
        const rbmd::Id pice_num,
        const rbmd::Id* atoms_type,
        const rbmd::Id* atoms_id,
        const rbmd::Id* start_id,
        const rbmd::Id* end_id,
        const rbmd::Id* id_verletlist,
        const rbmd::Id* id_random_neighbor,
        const rbmd::Id* random_neighbor_num,
        const Real7*  rhor_spline,
        const Real7*  z2r_spline,
        const rbmd::Real* px,
        const rbmd::Real* py,
        const rbmd::Real* pz,
        const  rbmd::Real* eam_fp,
        rbmd::Real* fx,
        rbmd::Real* fy,
        rbmd::Real* fz);
};


template <typename DEVICE>
struct ComputeEAMEnergy
{
  void operator()( Box box,
        EAMParameters eam_paras,
        const rbmd::Real cut_off,
        const rbmd::Id num_atoms,
        const rbmd::Id* atoms_type,
        const rbmd::Id* atoms_id,
        const rbmd::Id* start_id,
        const rbmd::Id* end_id,
        const rbmd::Id* id_verletlist,
        const Real7* rhor_spline,
        const Real7* frho_spline,
        const Real7* z2r_spline,
        const rbmd::Real* px,
        const rbmd::Real* py,
        const rbmd::Real* pz,
        rbmd::Real* eam_rho,
        rbmd::Real* eam_fp,
        rbmd::Real* energy_embedding,
        rbmd::Real* energy_pair);
};


//////////////////////////////////////////////
///

template <>
struct ComputeFpVerlet<device::DEVICE_GPU>
{
  void operator()( Box box,
        EAMParameters eam_paras,
        const rbmd::Real cut_off,
        const rbmd::Id num_atoms,
        const rbmd::Id* atoms_type,
        const rbmd::Id* atoms_id,
        const rbmd::Id* start_id,
        const rbmd::Id* end_id,
        const rbmd::Id* id_verletlist,
        const Real7*  rhor_spline,
        const Real7* frho_spline,
        const rbmd::Real* px,
        const rbmd::Real* py,
        const rbmd::Real* pz,
        rbmd::Real* eam_fp,
        rbmd::Real* energy_embedding);
};

template <>
struct ComputeEAMForceVerlet<device::DEVICE_GPU> {
  void operator()(Box box, EAMParameters eam_paras, const rbmd::Real cut_off,
                  const rbmd::Id num_atoms, const rbmd::Id* atoms_type,
                  const rbmd::Id* atoms_id, const rbmd::Id* start_id,
                  const rbmd::Id* end_id, const rbmd::Id* id_verletlist,
                  const Real7* rhor_spline, const Real7* frho_spline,
                  const Real7* z2r_spline, const rbmd::Real* px,
                  const rbmd::Real* py, const rbmd::Real* pz,
                  rbmd::Real* eam_rho, rbmd::Real* eam_fp, rbmd::Real* fx,
                  rbmd::Real* fy, rbmd::Real* fz, rbmd::Real* energy_embedding,
                  rbmd::Real* energy_pair);
};


template <>
struct ComputeRhoRBL<device::DEVICE_GPU>
{
  void operator()( Box box,
        EAMParameters eam_paras,
        const rbmd::Real rs,
        const rbmd::Real rc,
        const rbmd::Id num_atoms,
        const rbmd::Id neighbor_sample_num,
        const rbmd::Id pice_num,
        const rbmd::Id* atoms_type,
        const rbmd::Id* atoms_id,
        const rbmd::Id* start_id,
        const rbmd::Id* end_id,
        const rbmd::Id* id_verletlist,
        const rbmd::Id* id_random_neighbor,
        const rbmd::Id* random_neighbor_num,
        const Real7*  rhor_spline,
        const rbmd::Real* px,
        const rbmd::Real* py,
        const rbmd::Real* pz,
        rbmd::Real* eam_rho);
};



template <>
struct ComputeFpRBL<device::DEVICE_GPU>
{
  void operator()( Box box,
        EAMParameters eam_paras,
        const rbmd::Real rs,
        const rbmd::Real rc,
        const rbmd::Id num_atoms,
        const rbmd::Id neighbor_sample_num,
        const rbmd::Id pice_num,
        const rbmd::Id* atoms_type,
        const rbmd::Id* atoms_id,
        const rbmd::Id* start_id,
        const rbmd::Id* end_id,
        const rbmd::Id* id_verletlist,
        const rbmd::Id* id_random_neighbor,
        const rbmd::Id* random_neighbor_num,
        const Real7*  rhor_spline,
        const Real7* frho_spline,
        const rbmd::Real* px,
        const rbmd::Real* py,
        const rbmd::Real* pz,
        rbmd::Real* eam_fp);
};

template <>
struct ComputeEAMForceRBL<device::DEVICE_GPU>
{
  void operator()( Box box,
        EAMParameters eam_paras,
        const rbmd::Real rs,
        const rbmd::Real rc,
        const rbmd::Id num_atoms,
        const rbmd::Id neighbor_sample_num,
        const rbmd::Id pice_num,
        const rbmd::Id* atoms_type,
        const rbmd::Id* atoms_id,
        const rbmd::Id* start_id,
        const rbmd::Id* end_id,
        const rbmd::Id* id_verletlist,
        const rbmd::Id* id_random_neighbor,
        const rbmd::Id* random_neighbor_num,
        const Real7*  rhor_spline,
        const Real7*  z2r_spline,
        const rbmd::Real* px,
        const rbmd::Real* py,
        const rbmd::Real* pz,
        const  rbmd::Real* eam_fp,
        rbmd::Real* fx,
        rbmd::Real* fy,
        rbmd::Real* fz);
};

template <>
struct ComputeEAMEnergy<device::DEVICE_GPU>
{
  void operator()( Box box,
        EAMParameters eam_paras,
        const rbmd::Real cut_off,
        const rbmd::Id num_atoms,
        const rbmd::Id* atoms_type,
        const rbmd::Id* atoms_id,
        const rbmd::Id* start_id,
        const rbmd::Id* end_id,
        const rbmd::Id* id_verletlist,
        const Real7* rhor_spline,
        const Real7* frho_spline,
        const Real7* z2r_spline,
        const rbmd::Real* px,
        const rbmd::Real* py,
        const rbmd::Real* pz,
        rbmd::Real* eam_rho,
        rbmd::Real* eam_fp,
        rbmd::Real* energy_embedding,
        rbmd::Real* energy_pair);
};

}