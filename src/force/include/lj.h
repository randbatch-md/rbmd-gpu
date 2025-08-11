#pragma once
#include "../../common/types.h"
#include "force.h"
#include "model/box.h"
#include "neighbor_list/include/neighbor_list/neighbor_list.h"
#include "neighbor_list/include/neighbor_list_builder/full_neighbor_list_builder.h"
class LJ : public Force {
 public:
  LJ();
  virtual ~LJ();

 void Init() override;
 void  Execute() override;
 void EvaluatePotentialEnergy() override;

 void ComputeLJVerlet();
 void ComputeLJRBL() ;
 void ComputeLJEnergy();

 private:
  std::shared_ptr<BaseNeighborListBuilder> _rbl_neighbor_list_builder;
  std::shared_ptr<BaseNeighborListBuilder> _neighbor_list_builder;
  std::shared_ptr<NeighborList> _rbl_list;
  std::shared_ptr<NeighborList> _list;

  //energy
  rbmd::Real _e_vdwl= 0.0;
  rbmd::Real _e_pe = 0;
  rbmd::Real _e_pe_init = 0;

  rbmd::Real _e_vdwl_rbl= 0.0;
  rbmd::Real _e_pe_rbl = 0;

  rbmd::Real virial[6];

  //RBL
  std::string  _energy_rbl_flag = "yes"; //default
  std::string _neighbor_type;
  rbmd::Real _cut_off;

  rbmd::Real _corr_value_x =0.0;
  rbmd::Real _corr_value_y =0.0;
  rbmd::Real _corr_value_z =0.0;


};