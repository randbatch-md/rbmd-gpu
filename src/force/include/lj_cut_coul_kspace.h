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

  void ComputeQsqSum();
  rbmd::Real ComputeRMS(rbmd::Id kmax,rbmd::Real box_length,rbmd::Real q2);
  void SetKspacePara();
  void eik_dot_r();
  void coeffs();
  void  ComputeChargeStructureFactorEwald(
          Box box,
          rbmd::Id num_atoms,
          rbmd::Id Kmax,
          rbmd::Real alpha,
          rbmd::Real qqr2e,
          rbmd::Real* value_Re_array,
          rbmd::Real* value_Im_array);
  void ComputeEwlad();//Ewald

  void RBEInit(Box box,rbmd::Real alpha,rbmd::Id RBE_P);
  void GetPsampleKey();
  void ComputeChargeStructureFactorRBE(
         Box box,
         rbmd::Id num_atoms,
         rbmd::Id Kmax,
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
        rbmd::Id Kmax,
        rbmd::Real alpha,
        rbmd::Real qqr2e,
        rbmd::Real& ave_ekspace);   //Ewald  Energy

private:
  std::shared_ptr<BaseNeighborListBuilder> _rbl_neighbor_list_builder;
  std::shared_ptr<BaseNeighborListBuilder> _neighbor_list_builder;
  std::shared_ptr<NeighborList> _rbl_list;
  std::shared_ptr<NeighborList> _list;

  //energy
  rbmd::Real _ave_evdwl= 0.0;
  rbmd::Real _ave_ecoul= 0.0;
  rbmd::Real _ave_self_energy= 0.0;
  rbmd::Real _ave_ekspace= 0.0;
  rbmd::Real _ave_pe = 0;

  rbmd::Real _ave_evdwl_rbl = 0;
  rbmd::Real _ave_ecoul_rbl = 0;
  rbmd::Real _ave_pe_rbl = 0;

  //RBL
  std::string _neighbor_type;
  rbmd::Real _cut_off;

  rbmd::Real _corr_value_x = 0.0;
  rbmd::Real _corr_value_y = 0.0;
  rbmd::Real _corr_value_z = 0.0;

  //EWALD
  rbmd::Real _accuracy;
  rbmd::Real _g_ewald;
  rbmd::Real _alpha;
  rbmd::Real _q2;

  rbmd::Id _Kmax;
  rbmd::Id _Kmax3D;
  rbmd::Id kmax_x, kmax_y, kmax_z;
  Int3 _kmax_array;
  Real3 _unitk;
  rbmd::Real _gsqmx;
  rbmd::Id _kmax_x_orig, _kmax_y_orig, _kmax_z_orig;
  rbmd::Id kcount;

  std::vector<rbmd::Id> kxvecs,kyvecs,kzvecs; //_Kmax3D
  std::vector<rbmd::Real> ug;  //_Kmax3D
  std::vector<std::vector<rbmd::Real>> eg;  //_Kmax3D 3
  std::vector<std::vector<rbmd::Real>> vg; //_Kmax3D 6

  // std::vector<rbmd::Real> sfacrl,sfacrl_all;  //_Kmax3D
  // std::vector<rbmd::Real> sfacim,sfacim_all;  //_Kmax3D

  // std::vector<std::vector<std::vector<rbmd::Real>>> cs, sn;

  thrust::device_vector<rbmd::Real>  _d_cs,_d_sn;
  thrust::device_vector<rbmd::Real>  _d_sfacrl,_d_sfacrl_all;
  thrust::device_vector<rbmd::Real>  _d_sfacim,_d_sfacim_all;


  rbmd::Id _num_k;
  rbmd::Real* _h_Re_array;
  rbmd::Real* _h_Im_array;

  //RBE
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

