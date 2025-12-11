#pragma once
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include <memory>
#include <random>
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
  void ComputeChargeStructureFactorRBE_opt(
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

  //RBSOG
  void ComputeRBSOG();
  void RBSOGInit();
  void RBSOGSetup();
  void RBSOGSampleKSpace();
  void ComputeRBSOGFactor();

  //RBSOG functions
  rbmd::Real G_sigma(rbmd::Real sigma, rbmd::Real r) const;
  rbmd::Real Compute_W0(rbmd::Real r0, rbmd::Real b) const;
  rbmd::Real Gaussian(int kx, int ky, int kz, const Box& box, rbmd::Real sigma, rbmd::Real b, rbmd::Real w0, int Mmax, const thrust::host_vector<rbmd::Real>& coef) const;
  rbmd::Real Gaussian_modify(int kx, int ky, int kz, const Box& box, rbmd::Real sigma, rbmd::Real b, rbmd::Real w0, int Mmax, const thrust::host_vector<rbmd::Real>& coef_npt) const;
  rbmd::Real Gaussian_Fourier_Plus(rbmd::Real Kx, rbmd::Real Ky, rbmd::Real Kz, rbmd::Real sigma, rbmd::Real b, rbmd::Real w0, int Mmax, const thrust::host_vector<rbmd::Real>& coef) const;
  rbmd::Real Gaussian_Fourier_Plus_modify(rbmd::Real Kx, rbmd::Real Ky, rbmd::Real Kz, rbmd::Real sigma, rbmd::Real b, rbmd::Real w0, int Mmax, const thrust::host_vector<rbmd::Real>& coef_npt) const;
  rbmd::Real randn_box_muller( rbmd::Real Mean,  rbmd::Real SquareMargin);
  rbmd::Real MH_D_Modify(int xx, rbmd::Real factor) const;

  //
  std::string _coulomb_type;
  rbmd::Real  _cut_off;

  //EWALD
  double _accuracy;
  double _g_ewald;
  rbmd::Real _alpha;
  double _q2;
  double _sum_sq_charge;
  double  _sum_charge;
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

  //RBSOG - New members
  rbmd::Real _rbsog_b;
  rbmd::Real _rbsog_sigma;
  rbmd::Real _rbsog_omega;
  rbmd::Id   _rbsog_Mmax;
  rbmd::Id   _rbsog_Kcut;
  rbmd::Real _rbsog_w0;
  rbmd::Real _rbsog_S;
  rbmd::Real _rbsog_S_npt;

  // State saved at sampling time
  rbmd::Real _rbsog_S0;
  rbmd::Real _rbsog_S_npt0;
  Box        _rbsog_box0; // Box dimensions at sampling time

  // RBSOG coefficients
  thrust::host_vector<rbmd::Real>   _h_rbsog_sl;
  thrust::host_vector<rbmd::Real>   _h_rbsog_coef;
  thrust::host_vector<rbmd::Real>   _h_rbsog_coef_npt;
  thrust::host_vector<rbmd::Id> _idx_npt ; //

  // Device versions for kernels
  thrust::device_vector<rbmd::Real> _d_rbsog_sl;
  thrust::device_vector<rbmd::Real> _d_rbsog_coef;
  thrust::device_vector<rbmd::Real> _d_rbsog_coef_npt;

  // MCMC sampling state (host)
  thrust::host_vector<Int3>         _h_rbsog_K_Sample_int;
  thrust::host_vector<rbmd::Id>     _h_rbsog_K_Sample_x;
  thrust::host_vector<rbmd::Id>     _h_rbsog_K_Sample_y;
  thrust::host_vector<rbmd::Id>     _h_rbsog_K_Sample_z;
  thrust::host_vector<rbmd::Id>     _h_rbsog_idx_npt;

  // Sampled K-vectors (device)
  thrust::device_vector<rbmd::Real> _d_rbsog_K_Sample_x;
  thrust::device_vector<rbmd::Real> _d_rbsog_K_Sample_y;
  thrust::device_vector<rbmd::Real> _d_rbsog_K_Sample_z;
  thrust::device_vector<rbmd::Real> _d_rbsog_K_npt_x;
  thrust::device_vector<rbmd::Real> _d_rbsog_K_npt_y;
  thrust::device_vector<rbmd::Real> _d_rbsog_K_npt_z;
  thrust::device_vector<rbmd::Id>   _d_rbsog_idx_npt_all;

  rbmd::Real  xprd0,yprd0,zprd0;

  thrust::device_vector<rbmd::Real> _d_fac;
  thrust::device_vector<rbmd::Real> _d_fac_npt;

  // Random number generator for sampling
  std::mt19937 _gen;
};