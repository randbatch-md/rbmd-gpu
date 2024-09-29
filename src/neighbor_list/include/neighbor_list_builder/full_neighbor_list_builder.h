#pragma once
#include "base_neighbor_list_builder.h"
#include "../neighbor_list/neighbor_list.h"

class FullNeighborListBuilder : public BaseNeighborListBuilder {
public:
  explicit FullNeighborListBuilder();

  std::shared_ptr<NeighborList> Build() override;

protected:
  void ComputeNeighborCells() override;

  void ComputeNeighborCellsWithoutPBC() override;

  void EstimateNeighborsList() override;

  bool GenerateNeighborsList() override;

  std::shared_ptr<DeviceData> _device_data;
};



