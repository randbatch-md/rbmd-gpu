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
  thrust::device_vector<rbmd::Real> _d_total_evdwl;
  thrust::device_vector<rbmd::Real> _d_total_ecoul;

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
  rbmd::Real* _h_Re_array;
  rbmd::Real* _h_Im_array;

  //RBE
  std::string _coulomb_type;
  rbmd::Id _RBE_P;
  rbmd::Real _alpha;
  rbmd::Id _Kmax;
  rbmd::Real _qqr2e;
  rbmd::Id _num_k;
  thrust::device_vector<rbmd::Real>  _P_Sample_x;
  thrust::device_vector<rbmd::Real>  _P_Sample_y;
  thrust::device_vector<rbmd::Real>  _P_Sample_z;
  //
  thrust::device_vector<rbmd::Id>  _psample_key;
  thrust::device_vector<rbmd::Real> _rhok_real_redue;
  thrust::device_vector<rbmd::Real> _rhok_image_redue;

};

