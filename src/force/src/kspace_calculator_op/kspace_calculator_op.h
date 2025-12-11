#pragma once
#include "../../common/device_types.h"
#include "../../common/types.h"
#include "../../data_manager/include/model/box.h"
#include "../common/erf_table.h"

namespace op {

  template<typename DEVICE>
  struct GenerateIndexArrayOp
  {
    void operator()(
    const rbmd::Id  num_atoms,
    const rbmd::Id  RBE_P,
    rbmd::Id* psample_key);
  };

  template<typename DEVICE>
  struct SqchargeOp
  {
    void operator()(
      const rbmd::Id num_atoms,
      const rbmd::Real* charge,
      rbmd::Real* sq_charge);
  };

  template<typename DEVICE>
  struct SumchargeOp
  {
    void operator()(
      const rbmd::Id num_atoms,
      const rbmd::Real* charge,
      rbmd::Real* sum_sq_charge,
      rbmd::Real* sum_charge);
  };

  template <typename DEVICE>
  struct ComputeChargeStructureFactorOp
  {
    void operator()(
            const rbmd::Id num_atoms,
            const Real3 K,
            const rbmd::Real* charge,
            const rbmd::Real* px,
            const rbmd::Real* py,
            const rbmd::Real* pz,
            rbmd::Real* density_real,
            rbmd::Real* density_imag);
  };

  template <typename DEVICE>
  struct ComputeEwaldForceOp
  {
    void operator()(
             Box box,
            const rbmd::Id num_atoms,
            const Int3  Kmax,
            const rbmd::Real alpha,
            const rbmd::Real qqr2e,
            const rbmd::Real* real_array,
            const rbmd::Real* imag_array,
            const rbmd::Real* charge,
            const rbmd::Real* px,
            const rbmd::Real* py,
            const rbmd::Real* pz,
            rbmd::Real* fx,
            rbmd::Real* fy,
            rbmd::Real* fz,
            rbmd::Real* flat_virial);
  };

  template <typename DEVICE>
  struct ComputeRBEForceVirialOp
  {
    void operator()(
      Box box,const rbmd::Id P,
      const rbmd::Real alpha,const rbmd::Real qqrd2e,
      const rbmd::Real qsqsum,const rbmd::Real qsum,
      const rbmd::Real* real_array,const rbmd::Real* imag_array,
      const rbmd::Real* p_sample_x,const rbmd::Real* p_sample_y,
      const rbmd::Real* p_sample_z,rbmd::Real* virial_tensor,
      rbmd::Real* energy);
  };


  //RBE
  template <typename DEVICE>
  struct ComputePnumberChargeStructureFactorOp
  {
    void operator()(
        Box box,
       const rbmd::Id num_atoms,
       const rbmd::Id p_number,
       const rbmd::Real* charge,
       const rbmd::Real* p_sample_x,
       const rbmd::Real* p_sample_y,
       const rbmd::Real* p_sample_z,
       const rbmd::Real* px,
       const rbmd::Real* py,
       const rbmd::Real* pz,
       rbmd::Real* density_real,
       rbmd::Real* density_imag);
  };

  template <typename DEVICE>
  struct ComputeRBEForceOp
  {
    void operator()(
      Box box,const rbmd::Id num_atoms,
     const rbmd::Id  p_number,const rbmd::Real alpha,
     const rbmd::Real qqr2e,const rbmd::Real* real_array,
     const rbmd::Real* imag_array,const rbmd::Real* charge,
     const rbmd::Real* p_sample_x,const rbmd::Real* p_sample_y,
     const rbmd::Real* p_sample_z,const rbmd::Real* px,
     const rbmd::Real* py,const rbmd::Real* pz,
     rbmd::Real* fx,rbmd::Real* fy,rbmd::Real* fz);
  };


  template <typename DEVICE>
  struct ComputeRBSOGFactorOp
  {
    void operator()(
    Box box, const rbmd::Id P, const rbmd::Real sigma, const rbmd::Real b, const rbmd::Id Mmax,
    const rbmd::Real L_ratio, const rbmd::Real S_ratio,
    const rbmd::Real S_npt_ratio,
    const rbmd::Real* K_x, const rbmd::Real* K_y, const rbmd::Real* K_z,
    const rbmd::Real* K_npt_x, const rbmd::Real* K_npt_y, const rbmd::Real* K_npt_z,
    const rbmd::Real* coef, const rbmd::Real* coef_npt,
    rbmd::Real* fac, rbmd::Real* fac_npt);
  };

  template <typename DEVICE>
  struct ComputePnumberChargeStructureFactorSOGOp
  {
    void operator()(
      Box box,const rbmd::Id num_atoms,const rbmd::Id p_number,
     const rbmd::Real* charge,const rbmd::Real* p_sample_x,
     const rbmd::Real* p_sample_y,const rbmd::Real* p_sample_z,
     const rbmd::Real* px,const rbmd::Real* py,const rbmd::Real* pz,
     rbmd::Real* density_real,rbmd::Real* density_imag);
  };

  template <typename DEVICE>
  struct ComputeRBSOGSampleForceOp
  {
    void operator()(
    Box box, const rbmd::Id num_atoms, const rbmd::Id P,
    const rbmd::Real qqr2e, const rbmd::Real S0_sample, const rbmd::Real S_npt_sample,
    const rbmd::Real* K_x, const rbmd::Real* K_y, const rbmd::Real* K_z,
    const rbmd::Real* K_npt_x, const rbmd::Real* K_npt_y, const rbmd::Real* K_npt_z,
    const rbmd::Id* idx_npt_all,const rbmd::Real* fac, const rbmd::Real* fac_npt,
    const rbmd::Real* density_real, const rbmd::Real* density_imag,
    const rbmd::Real* charge,const rbmd::Real* px, const rbmd::Real* py, const rbmd::Real* pz,
    rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,rbmd::Real* global_virial,
    rbmd::Real* energy_parts);
  };

  template <typename DEVICE>
  struct ComputeDirectChargeStructureFactorOp
  {
    void operator()(
    const rbmd::Id num_atoms, const rbmd::Id num_k_direct,
    const rbmd::Real*  k_direct_x,  const rbmd::Real*  k_direct_y,
    const rbmd::Real*  k_direct_z,const rbmd::Real* charge,
    const rbmd::Real* px,const rbmd::Real* py,
    const rbmd::Real* pz,rbmd::Real* density_real,rbmd::Real* density_imag);
  };

  template <typename DEVICE>
  struct ComputeRBSOGDirectForceOp
  {
    void operator()(
    Box box, const rbmd::Id num_atoms, const rbmd::Id num_k_direct,
    const rbmd::Real qqr2e,    const rbmd::Real*  k_direct_x,
    const rbmd::Real*  k_direct_y,const rbmd::Real*  k_direct_z,
    const rbmd::Real* f_b_sigma,const rbmd::Real* f_b_sigma_npt,
    const rbmd::Real* density_real,const rbmd::Real* density_imag,
    const rbmd::Real* charge,const rbmd::Real* px, const rbmd::Real* py, const rbmd::Real* pz,
    rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,rbmd::Real* global_virial,
    rbmd::Real* energy_parts);
  };


/////////////////////////////////////////

template<>
struct GenerateIndexArrayOp<device::DEVICE_GPU>
{
  void operator()(
  const rbmd::Id  num_atoms,
  const rbmd::Id  RBE_P,
  rbmd::Id* psample_key);
};

template<>
struct SqchargeOp<device::DEVICE_GPU>
{
  void operator()(
    const rbmd::Id num_atoms,
    const rbmd::Real* charge,
    rbmd::Real* sq_charge);
};

template<>
struct SumchargeOp<device::DEVICE_GPU>
{
  void operator()(
    const rbmd::Id num_atoms,
    const rbmd::Real* charge,
    rbmd::Real* sum_sq_charge,
    rbmd::Real* sum_charge);
};

template <>
struct ComputeChargeStructureFactorOp<device::DEVICE_GPU>
{
  void operator()(
          const rbmd::Id num_atoms,
          const Real3 K,
          const rbmd::Real* charge,
          const rbmd::Real* px,
          const rbmd::Real* py,
          const rbmd::Real* pz,
          rbmd::Real* density_real,
          rbmd::Real* density_imag);
};

template <>
struct ComputeEwaldForceOp<device::DEVICE_GPU>
{
  void operator()(
           Box box,
          const rbmd::Id num_atoms,
          const Int3  Kmax,
          const rbmd::Real alpha,
          const rbmd::Real qqr2e,
          const rbmd::Real* real_array,
          const rbmd::Real* imag_array,
          const rbmd::Real* charge,
          const rbmd::Real* px,
          const rbmd::Real* py,
          const rbmd::Real* pz,
          rbmd::Real* fx,
          rbmd::Real* fy,
          rbmd::Real* fz,
          rbmd::Real* flat_virial);
};

template <>
struct ComputePnumberChargeStructureFactorOp<device::DEVICE_GPU>
{
  void operator()(
      Box box,
     const rbmd::Id num_atoms,
     const rbmd::Id p_number,
     const rbmd::Real* charge,
     const rbmd::Real* p_sample_x,
     const rbmd::Real* p_sample_y,
     const rbmd::Real* p_sample_z,
     const rbmd::Real* px,
     const rbmd::Real* py,
     const rbmd::Real* pz,
     rbmd::Real* density_real,
     rbmd::Real* density_imag);
};

template <>
struct ComputeRBEForceOp<device::DEVICE_GPU>
{
  void operator()(
      Box box,const rbmd::Id num_atoms,
     const rbmd::Id  p_number,const rbmd::Real alpha,
     const rbmd::Real qqr2e,const rbmd::Real* real_array,
     const rbmd::Real* imag_array,const rbmd::Real* charge,
     const rbmd::Real* p_sample_x,const rbmd::Real* p_sample_y,
     const rbmd::Real* p_sample_z,const rbmd::Real* px,
     const rbmd::Real* py,const rbmd::Real* pz,
     rbmd::Real* fx,rbmd::Real* fy,rbmd::Real* fz);
};

  template <>
  struct ComputeRBEForceVirialOp<device::DEVICE_GPU>
  {
    void operator()(
    Box box,const rbmd::Id P,
    const rbmd::Real alpha,const rbmd::Real qqrd2e,
    const rbmd::Real qsqsum,const rbmd::Real qsum,
    const rbmd::Real* real_array,const rbmd::Real* imag_array,
    const rbmd::Real* p_sample_x,const rbmd::Real* p_sample_y,
    const rbmd::Real* p_sample_z,rbmd::Real* virial_tensor,
    rbmd::Real* energy);
  };

  template <>
  struct ComputeRBSOGFactorOp<device::DEVICE_GPU>
  {
    void operator()(
    Box box, const rbmd::Id P, const rbmd::Real sigma, const rbmd::Real b, const rbmd::Id Mmax,
    const rbmd::Real L_ratio, const rbmd::Real S_ratio,
    const rbmd::Real S_npt_ratio,
    const rbmd::Real* K_x, const rbmd::Real* K_y, const rbmd::Real* K_z,
    const rbmd::Real* K_npt_x, const rbmd::Real* K_npt_y, const rbmd::Real* K_npt_z,
    const rbmd::Real* coef, const rbmd::Real* coef_npt,
    rbmd::Real* fac, rbmd::Real* fac_npt);
  };

  template <>
  struct ComputePnumberChargeStructureFactorSOGOp<device::DEVICE_GPU>
  {
    void operator()(
        Box box,const rbmd::Id num_atoms,const rbmd::Id p_number,
       const rbmd::Real* charge,const rbmd::Real* p_sample_x,
       const rbmd::Real* p_sample_y,const rbmd::Real* p_sample_z,
       const rbmd::Real* px,const rbmd::Real* py,const rbmd::Real* pz,
       rbmd::Real* density_real,rbmd::Real* density_imag);
  };

  template <>
  struct ComputeRBSOGSampleForceOp<device::DEVICE_GPU>
  {
    void operator()(
    Box box, const rbmd::Id num_atoms, const rbmd::Id P,
    const rbmd::Real qqr2e, const rbmd::Real S0_sample, const rbmd::Real S_npt_sample,
    const rbmd::Real* K_x, const rbmd::Real* K_y, const rbmd::Real* K_z,
    const rbmd::Real* K_npt_x, const rbmd::Real* K_npt_y, const rbmd::Real* K_npt_z,
    const rbmd::Id* idx_npt_all,const rbmd::Real* fac, const rbmd::Real* fac_npt,
    const rbmd::Real* density_real, const rbmd::Real* density_imag,
    const rbmd::Real* charge,const rbmd::Real* px, const rbmd::Real* py, const rbmd::Real* pz,
    rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,rbmd::Real* global_virial,
    rbmd::Real* energy_parts);
  };

  template <>
  struct ComputeDirectChargeStructureFactorOp<device::DEVICE_GPU>
  {
    void operator()(
    const rbmd::Id num_atoms, const rbmd::Id num_k_direct,
    const rbmd::Real*  k_direct_x,  const rbmd::Real*  k_direct_y,
    const rbmd::Real*  k_direct_z,const rbmd::Real* charge,
    const rbmd::Real* px,const rbmd::Real* py,
    const rbmd::Real* pz,rbmd::Real* density_real,rbmd::Real* density_imag);
  };

  template <>
  struct ComputeRBSOGDirectForceOp<device::DEVICE_GPU>
  {
    void operator()(
    Box box, const rbmd::Id num_atoms, const rbmd::Id num_k_direct,
    const rbmd::Real qqr2e,    const rbmd::Real*  k_direct_x,
    const rbmd::Real*  k_direct_y,const rbmd::Real*  k_direct_z,
    const rbmd::Real* f_b_sigma,const rbmd::Real* f_b_sigma_npt,
    const rbmd::Real* density_real,const rbmd::Real* density_imag,
    const rbmd::Real* charge,const rbmd::Real* px, const rbmd::Real* py, const rbmd::Real* pz,
    rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz,rbmd::Real* global_virial,
    rbmd::Real* energy_parts);
  };

}