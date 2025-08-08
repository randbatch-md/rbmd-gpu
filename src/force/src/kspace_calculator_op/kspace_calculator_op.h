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
    const Box box,
    const rbmd::Id P,             // 随机采样波矢数量
    const rbmd::Real alpha,       // Ewald参数
    const rbmd::Real qqrd2e,
    const rbmd::Real qsqsum,
    const rbmd::Real qsum,
    const rbmd::Real* real_array, // ρ(k)实部数组
    const rbmd::Real* imag_array, // ρ(k)虚部数组
    const rbmd::Real* p_sample_x, // 采样波矢x分量
    const rbmd::Real* p_sample_y, // 采样波矢y分量
    const rbmd::Real* p_sample_z, // 采样波矢z分量
    rbmd::Real* virial_tensor,
    rbmd::Real* energy);    // 输出维里张量[6]);
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
        Box box,
       const rbmd::Id num_atoms,
       const rbmd::Id  p_number,
       const rbmd::Real alpha,
       const rbmd::Real qqr2e,
       const rbmd::Real* real_array,
       const rbmd::Real* imag_array,
       const rbmd::Real* charge,
       const rbmd::Real* p_sample_x,
       const rbmd::Real* p_sample_y,
       const rbmd::Real* p_sample_z,
       const rbmd::Real* px,
       const rbmd::Real* py,
       const rbmd::Real* pz,
       rbmd::Real* fx,
       rbmd::Real* fy,
       rbmd::Real* fz,
       rbmd::Real* flat_virial);
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
      Box box,
     const rbmd::Id num_atoms,
     const rbmd::Id  p_number,
     const rbmd::Real alpha,
     const rbmd::Real qqr2e,
     const rbmd::Real* real_array,
     const rbmd::Real* imag_array,
     const rbmd::Real* charge,
     const rbmd::Real* p_sample_x,
     const rbmd::Real* p_sample_y,
     const rbmd::Real* p_sample_z,
     const rbmd::Real* px,
     const rbmd::Real* py,
     const rbmd::Real* pz,
     rbmd::Real* fx,
     rbmd::Real* fy,
     rbmd::Real* fz,
     rbmd::Real* flat_virial);
};

  template <>
  struct ComputeRBEForceVirialOp<device::DEVICE_GPU>
  {
    void operator()(
    const Box box,
    const rbmd::Id P,             // 随机采样波矢数量
    const rbmd::Real alpha,       // Ewald参数
    const rbmd::Real qqrd2e,
    const rbmd::Real qsqsum,
    const rbmd::Real qsum,
    const rbmd::Real* real_array, // ρ(k)实部数组
    const rbmd::Real* imag_array, // ρ(k)虚部数组
    const rbmd::Real* p_sample_x, // 采样波矢x分量
    const rbmd::Real* p_sample_y, // 采样波矢y分量
    const rbmd::Real* p_sample_z, // 采样波矢z分量
    rbmd::Real* virial_tensor ,
    rbmd::Real* energy);    // 输出维里张量[6]);
  };


}