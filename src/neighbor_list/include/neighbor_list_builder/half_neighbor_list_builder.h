#pragma once
#include "full_neighbor_list_builder.h"

class HalfNeighborListBuilder : public FullNeighborListBuilder {
public:
  explicit HalfNeighborListBuilder();

protected:
  void ComputeNeighborCells() override;

  void EstimateNeighborsList() override;

  bool GenerateNeighborsList() override;

  bool _without_pbc = false;

};