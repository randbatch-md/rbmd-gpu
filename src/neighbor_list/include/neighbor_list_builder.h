#pragma once
#include "box.h"
#include "common/object.h"
#include "common/types.h"
#include "linked_cell.h"
#include "neighbor_list.h"

// linked Cell should new in out
class NeighborListBuilder : public Object {
public:
  explicit NeighborListBuilder(std::shared_ptr<LinkedCell> linked_cell);

  virtual std::shared_ptr<NeighborList> Build() = 0;

protected:
  std::shared_ptr<LinkedCell> _linked_cell =
      nullptr; // linkedcell应该很早就开始创建
  std::shared_ptr<NeighborList> _neighbor_list = nullptr;

  virtual void ComputeNeighborCells() = 0;

  virtual void EstimateNeighbousList() = 0;

  virtual bool GenerateNeighbousList() = 0;

  void ReductionSum(rbmd::Id* d_src_array, rbmd::Id* d_dst, rbmd::Id size);

  void InitNeighborListIndices();

  rbmd::Id _neighbor_cell_num = 0;
  bool should_realloc = true;
  Box* _d_box; // TODO get by devie_data
};