#pragma once
#include "neighbor_list_builder.h"

class FullNeighborListBuilder : public NeighborListBuilder {
public:
  explicit FullNeighborListBuilder(
      const std::shared_ptr<LinkedCell>& linked_cell);

  std::shared_ptr<NeighborList> Build() override;

protected:
  void ComputeNeighborCells() override;

  void EstimateNeighbousList() override;

  bool GenerateNeighbousList() override;

  std::shared_ptr<DeviceData> _device_data;
};



