#pragma once
#include "../../common/types.h"
#include "force.h"
#include "model/box.h"
#include "neighbor_list/include/neighbor_list/neighbor_list.h"
#include "neighbor_list/include/neighbor_list_builder/full_neighbor_list_builder.h"
class LJCoulEwaldForce : public Force
{
public:
  LJCoulEwaldForce();
  virtual ~LJCoulEwaldForce();

  void Init() override;
  void  Execute() override;

  void ComputeLJCutCoulForce();
  //
  void  ComputeChargeStructureFactorEwald(
          Box* box,
          rbmd::Id _num_atoms,
          rbmd::Id Kmax,
          rbmd::Real alpha,
          rbmd::Real qqr2e,
          rbmd::Real* value_Re_array,
          rbmd::Real* value_Im_array);

  void ComputeEwladForce();
  void RBEInit(rbmd::Real alpha,Box* box);
  void ComputeRBEForce();
  void SumForces();
  void ComputeSelfEnergy(
    rbmd::Real alpha,
    rbmd::Real qqr2e);

private:
  rbmd::Id _num_atoms;
  std::shared_ptr<BaseNeighborListBuilder> _rbl_neighbor_list_builder;
  std::shared_ptr<BaseNeighborListBuilder> _neighbor_list_builder;
  std::shared_ptr<NeighborList> rbl_list;
  std::shared_ptr<NeighborList> list;
  Box box;

  rbmd::Real _corr_value_x;
  rbmd::Real _corr_value_y;
  rbmd::Real _corr_value_z;
  rbmd::Real* _d_total_evdwl;
  rbmd::Real* _d_total_ecoul;


  //energy
  rbmd::Real _ave_self_energy;
  rbmd::Real _ave_eewald;
  rbmd::Id _num_k;

  //
  rbmd::Id _RBE_P;
  rbmd::Real* _P_Sample_x;
  rbmd::Real* _P_Sample_y;
  rbmd::Real* _P_Sample_z;

  //
  thrust::device_vector<rbmd::Id>  _psample_key;
  thrust::device_vector<rbmd::Real>  _rhok_real;
  thrust::device_vector<rbmd::Real>  _rhok_image;
};

