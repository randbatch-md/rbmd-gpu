#pragma once
#include "../../data_manager/include/model/box.h"
#include "../common/object.h"
#include "../common/types.h"
#include "../neighbor_list/include/linked_cell/linked_cell.h"
#include "../neighbor_list/include/neighbor_list.h"

// linked Cell should new in out
class NeighborListBuilder : public Object {
 public:
  explicit NeighborListBuilder();
  ~NeighborListBuilder() override;

  virtual std::shared_ptr<NeighborList> Build() = 0;

 protected:
  std::shared_ptr<LinkedCell> _linked_cell;
  std::shared_ptr<NeighborList> _neighbor_list = nullptr;

  virtual void ComputeNeighborCells() = 0;

  virtual void ComputeNeighborCellsWithoutPBC() = 0;

  virtual void EstimateNeighborsList() = 0;

  virtual bool GenerateNeighborsList() = 0;

  void ReductionSum(rbmd::Id* d_src_array, rbmd::Id* d_dst, rbmd::Id size);

  void InitNeighborListIndices();

  rbmd::Id _neighbor_cell_num = 0;
  bool should_realloc = true;
  Box* _d_box;   // TODO 可能不太适合 待重构
  bool* _d_should_realloc;

};
