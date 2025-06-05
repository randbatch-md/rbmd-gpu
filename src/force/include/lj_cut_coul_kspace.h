#pragma once
#include "../../common/types.h"
#include "force.h"
#include "model/box.h"
#include "../common/erf_table.h"
#include "neighbor_list/include/neighbor_list/neighbor_list.h"
#include "neighbor_list/include/neighbor_list_builder/full_neighbor_list_builder.h"
class LJCutCoulKspace : public Force
{
public:
  LJCutCoulKspace();
  virtual ~LJCutCoulKspace();

  void Init() override;
  void  Execute() override;
  void EvaluatePotentialenergy() override;

  void ComputeLJCutCoulForce();
  void ComputeLJVerlet();
  void ComputeLJRBL();

  void ComputeKspaceForce();
  void SumForces();

  //
  void ComputeQsqSum();
  rbmd::Real ComputeRMS(rbmd::Id kmax,rbmd::Real box_length,rbmd::Real q2);
  void SetKspacePara();
  void ComputeQsf();//Charge Structure Factor
  void ComputeQsf_fix();
  void coeffs();
  void ComputeEwlad_fix();

  void ComputeWaveVectors();
  void ComputeWaveVectors_2( Box box,rbmd::Id Kmax);

  void  ComputeChargeStructureFactorEwald(
          Box box,
          rbmd::Id num_atoms,
          Int3 Kmax_array,
          rbmd::Real alpha,
          rbmd::Real qqr2e,
          thrust::host_vector<rbmd::Real> value_Re_array,
          thrust::host_vector<rbmd::Real> value_Im_array);
  void ComputeEwlad();//Ewald

  void RBEInit(Box box,rbmd::Real alpha,rbmd::Id RBE_P);
  void GetPsampleKey();
  void ComputeChargeStructureFactorRBE(
         Box box,
         rbmd::Id num_atoms,
         Int3 Kmax_array,
         rbmd::Real alpha,
         rbmd::Id RBE_P,
         rbmd::Real qqr2e,
         thrust::device_vector<rbmd::Real> rhok_real_redue,
         thrust::device_vector<rbmd::Real> rhok_image_redue);
  void ComputeRBE();//RBE


  void ComputeLJCoulEnergy();

  void ComputeSelfEnergy(
    rbmd::Real alpha,
    rbmd::Real qqr2e,
    rbmd::Real& ave_self_energy);  //self  Energy

  void ComputeKspaceEnergy(
        Box box,
        rbmd::Id _num_atoms,
        Int3 Kmax_array,
        rbmd::Real alpha,
        rbmd::Real qqr2e,
        rbmd::Real& ave_ekspace);   //Ewald  Energy

private:
  std::shared_ptr<BaseNeighborListBuilder> _rbl_neighbor_list_builder;
  std::shared_ptr<BaseNeighborListBuilder> _neighbor_list_builder;
  std::shared_ptr<NeighborList> _rbl_list;
  std::shared_ptr<NeighborList> _list;

  //energy
  rbmd::Real _e_vdwl= 0.0;
  rbmd::Real _e_coul= 0.0;
  rbmd::Real _e_self_energy= 0.0;
  rbmd::Real _e_kspace= 0.0;
  rbmd::Real _e_pe = 0;

  rbmd::Real _e_vdwl_rbl = 0;
  rbmd::Real _e_coul_rbl = 0;
  rbmd::Real _e_pe_rbl = 0;

  //RBL
  std::string _neighbor_type;
  rbmd::Real _cut_off;
  std::string  _energy_rbl_flag = "yes";

  rbmd::Real _corr_value_x = 0.0;
  rbmd::Real _corr_value_y = 0.0;
  rbmd::Real _corr_value_z = 0.0;

  //EWALD
  rbmd::Real _accuracy;
  rbmd::Real _g_ewald;
  rbmd::Real _alpha;
  rbmd::Real _sum_sq_charge;

  rbmd::Id _Kmax;
  rbmd::Id _Kmax3D;
  rbmd::Id kmax_x, kmax_y, kmax_z;
  Int3 _kmax_array;
  Real3 _unitk;
  rbmd::Real _gsqmx;
  rbmd::Id _kmax_x_orig, _kmax_y_orig, _kmax_z_orig;
  rbmd::Id kcount;

  std::vector<rbmd::Id> kxvecs,kyvecs,kzvecs; //_Kmax3D
  std::vector<rbmd::Id> kxvecs_R,kyvecs_R,kzvecs_R; //_Kmax3D
  std::vector<rbmd::Real> ug;  //_Kmax3D
  std::vector<Int3> kmax_vec3D; //_Kmax3D 3

  // std::vector<rbmd::Real> eg_flat;
  // std::vector<rbmd::Real> vg_flat;
  thrust::host_vector<rbmd::Real> eg_flat;
  thrust::host_vector<rbmd::Real> vg_flat;

  thrust::device_vector<Int3> _d_kmax_vec3D; //_Kmax3D 3
  thrust::device_vector<rbmd::Real> _d_eg_flat;
  thrust::device_vector<rbmd::Real> _d_vg_flat;

  // std::vector<rbmd::Real> sfacrl,sfacrl_all;  //_Kmax3D
  // std::vector<rbmd::Real> sfacim,sfacim_all;  //_Kmax3D

  // std::vector<std::vector<std::vector<rbmd::Real>>> cs, sn;

  thrust::device_vector<rbmd::Real>  _d_cs,_d_sn;
  thrust::device_vector<rbmd::Real>  _d_qfactor_real,_d_qfactor_real_all;
  thrust::device_vector<rbmd::Real>  _d_qfactor_image,_d_qfactor_image_all;


  thrust::host_vector<rbmd::Real> wavevec_indices_host;
  thrust::device_vector<rbmd::Id>  _d_kxvecs,_d_kyvecs,_d_kzvecs;
  thrust::device_vector<rbmd::Real>  _d_kxvecs_R,_d_kyvecs_R,_d_kzvecs_R;

  rbmd::Id _num_k;
  rbmd::Real* _h_Re_array;
  rbmd::Real* _h_Im_array;

  //RBE
  std::string  _energy_rbe_flag = "yes";
  std::string _coulomb_type;
  rbmd::Id _RBE_P;
  rbmd::Real _qqr2e;

  thrust::device_vector<rbmd::Real>  _P_Sample_x;
  thrust::device_vector<rbmd::Real>  _P_Sample_y;
  thrust::device_vector<rbmd::Real>  _P_Sample_z;
  //
  thrust::device_vector<rbmd::Id>  _psample_key;
  thrust::device_vector<rbmd::Real> _rhok_real_redue;
  thrust::device_vector<rbmd::Real> _rhok_image_redue;

};

