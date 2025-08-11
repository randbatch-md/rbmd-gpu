#pragma once
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include <memory>
#include <string>

#include "../../common/types.h"
#include "force.h"

class KSpaceCalculator : public Force{
public:
  KSpaceCalculator();

  virtual ~KSpaceCalculator() = default;

  void Init() override;
  void Execute() override;
  void EvaluatePotentialEnergy() override ;

  rbmd::Real GetKspacEnergy() const { return _e_kspace; }

private:
  void ComputeQsqSum();
  void GetPsampleKey();
  rbmd::Real ComputeRMS(rbmd::Id kmax,rbmd::Real box_length,rbmd::Real q2);
  void SetKspacePara();

  //Ewald
  void  ComputeChargeStructureFactorEwald(
        Box box,
        rbmd::Id num_atoms,
        Int3 Kmax_array,
        rbmd::Real alpha,
        rbmd::Real qqr2e,
        thrust::host_vector<rbmd::Real>& value_Re_array,
        thrust::host_vector<rbmd::Real>& value_Im_array);
  void ComputeEwald();

  //RBE
  void RBEInit(Box box,rbmd::Real alpha,rbmd::Id RBE_P);
  void ComputeChargeStructureFactorRBE(
       Box box,
       rbmd::Id num_atoms,
       Int3 Kmax_array,
       rbmd::Real alpha,
       rbmd::Id RBE_P,
       rbmd::Real qqr2e,
       thrust::device_vector<rbmd::Real>& rhok_real_redue,
       thrust::device_vector<rbmd::Real>& rhok_image_redue);
  void ComputeRBE();
  void ComputeRBEVirial();

  void ComputeSelfEnergy(rbmd::Real alpha,rbmd::Real qqr2e,rbmd::Real& ave_self_energy);  //self  Energy
  void ComputeKspaceEnergy(
        Box box,
        rbmd::Id _num_atoms,
        Int3 Kmax_array,
        rbmd::Real alpha,
        rbmd::Real qqr2e,
        rbmd::Real& ave_ekspace);   //Ewald  Energy

  //
  std::string _coulomb_type;
  rbmd::Real  _cut_off;

  //EWALD
  rbmd::Real _accuracy;
  rbmd::Real _g_ewald;
  rbmd::Real _alpha;
  rbmd::Real _q2;
  rbmd::Real _sum_sq_charge;
  rbmd::Real  _sum_charge;
  rbmd::Real* _h_Re_array;
  rbmd::Real* _h_Im_array;

  rbmd::Id _Kmax;
  rbmd::Id _Kmax3D;
  rbmd::Id kmax_x, kmax_y, kmax_z;
  Int3 _kmax_array;
  Real3 _unitk;
  rbmd::Real _gsqmx;
  rbmd::Id _kmax_x_orig, _kmax_y_orig, _kmax_z_orig;
  rbmd::Id kcount;

  //RBE
  rbmd::Id _RBE_P;
  std::string  _energy_rbe_flag = "yes"; //default
  rbmd::Real _qqr2e;
  rbmd::Id _num_k;
  thrust::device_vector<rbmd::Real>  _P_Sample_x;
  thrust::device_vector<rbmd::Real>  _P_Sample_y;
  thrust::device_vector<rbmd::Real>  _P_Sample_z;
  //
  thrust::device_vector<rbmd::Id>  _psample_key;
  thrust::device_vector<rbmd::Real> _rhok_real_redue;
  thrust::device_vector<rbmd::Real> _rhok_image_redue;
  //
  rbmd::Real _e_kspace = 0;
  rbmd::Real _e_self_energy = 0;

};