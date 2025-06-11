#pragma once
#include "../neighbor_list/neighbor_list.h"
#include "base_neighbor_list_builder.h"

class MACENeighborListBuilder : public BaseNeighborListBuilder {
 public:
  explicit MACENeighborListBuilder();

  std::shared_ptr<NeighborList> Build() override;

 protected:
  void ComputeNeighborCells() override;

  void ComputeNeighborCellsWithoutPBC() override;

  void EstimateNeighborsList() override;

  rbmd::Id GenerateNeighborsList() override;

  std::shared_ptr<DeviceData> _device_data;
};
