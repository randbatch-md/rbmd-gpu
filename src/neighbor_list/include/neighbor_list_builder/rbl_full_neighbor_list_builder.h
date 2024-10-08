#pragma once
#include "data_manager.h"
#include "full_neighbor_list_builder.h"

class RblFullNeighborListBuilder :public FullNeighborListBuilder {
public:
  explicit RblFullNeighborListBuilder();
protected:
  rbmd::Real _r_core = 0;
  rbmd::Id _neighbor_sample_num = 0;
  rbmd::Real _system_rho = 0;
  rbmd::Id _selection_frequency = 0;

  rbmd::Id GenerateNeighborsList() override;
  void GetRblParams();


};