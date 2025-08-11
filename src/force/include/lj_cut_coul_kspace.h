#pragma once
#include "../../common/types.h"
#include "../common/erf_table.h"
#include "force.h"
#include "kspace_calculator.h"
#include "model/box.h"
#include "neighbor_list/include/neighbor_list/neighbor_list.h"
#include "neighbor_list/include/neighbor_list_builder/full_neighbor_list_builder.h"
class LJCutCoulKspace : public Force
{
public:
  LJCutCoulKspace();
  virtual ~LJCutCoulKspace();

  void Init() override;
  void  Execute() override;
  void EvaluatePotentialEnergy() override;

  void ComputeLJCutCoulForce();
  void ComputeLJVerlet();
  void ComputeLJRBL();
  void ComputeLJCoulEnergy();

  void SumForces();

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
  std::string  _energy_rbl_flag = "yes"; //default

  rbmd::Real _corr_value_x = 0.0;
  rbmd::Real _corr_value_y = 0.0;
  rbmd::Real _corr_value_z = 0.0;

  //kspace
  rbmd::Real _qqr2e;
  rbmd::Real _alpha;
  std::unique_ptr<KSpaceCalculator> _kspace_calculator;
};

