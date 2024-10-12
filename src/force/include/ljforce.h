#pragma once
#include "../../common/types.h"
#include "force.h"
#include "model/box.h"
#include "neighbor_list/include/neighbor_list/neighbor_list.h"
#include "neighbor_list/include/neighbor_list_builder/full_neighbor_list_builder.h"
class LJForce : public Force {
 public:
  LJForce();
  virtual ~LJForce();

	void Init() override;
	void  Execute() override;

 private:
  rbmd::Id _num_atoms;
  std::shared_ptr<BaseNeighborListBuilder> _rbl_neighbor_list_builder;
  std::shared_ptr<BaseNeighborListBuilder> _neighbor_list_builder;
  std::shared_ptr<NeighborList> _rbl_list;
  std::shared_ptr<NeighborList> _list;
  Box box;

  rbmd::Real _corr_value_x;
  rbmd::Real _corr_value_y;
  rbmd::Real _corr_value_z;
  rbmd::Real* _d_total_evdwl;

  //RBL
  std::string _neighbor_type;
  rbmd::Real _cut_off;

  rbmd::Real virial[6];
};