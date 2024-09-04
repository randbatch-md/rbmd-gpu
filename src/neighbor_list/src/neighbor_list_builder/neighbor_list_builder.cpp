#include "neighbor_list_builder.h"

#include <hipcub/hipcub.hpp>

#include "common/device_types.h"
#include "neighbor_list/src/linked_cell_op/neighbor_list_op.h"

NeighborListBuilder::NeighborListBuilder(
    std::shared_ptr<LinkedCell> linked_cell) {
  this->_linked_cell = linked_cell;
  // todo neighbor num by Cell in cutoff
}

void NeighborListBuilder::ReductionSum(rbmd::Id* d_src_array, rbmd::Id* d_dst,
                                       rbmd::Id size) {
  void* temp = nullptr;
  size_t temp_bytes = 0;
  CHECK_RUNTIME(hipcub::DeviceReduce::Sum(temp, temp_bytes, d_src_array, d_dst,
    static_cast<int>(size)));
  CHECK_RUNTIME(MALLOC(&temp, temp_bytes));
  CHECK_RUNTIME(hipcub::DeviceReduce::Sum(temp, temp_bytes, d_src_array, d_dst,
    static_cast<int>(size)));
  CHECK_RUNTIME(FREE(temp));
}

void NeighborListBuilder::InitNeighborListIndices() {
  void* temp = nullptr;
  size_t temp_bytes = 0;

  CHECK_RUNTIME(hipcub::DeviceScan::ExclusiveSum(
    temp, temp_bytes,
    thrust::raw_pointer_cast(_neighbor_list->_d_max_neighbor_num.data()),
    thrust::raw_pointer_cast(_neighbor_list->_start_idx.data()),
    static_cast<int>(_linked_cell->_total_atoms_num)));
  CHECK_RUNTIME(MALLOC(&temp, temp_bytes));
  // 重新设置开始索引   start按照max为间隔，后面加上前面的
  // https://blog.csdn.net/qq_45914558/article/details/107385862
  CHECK_RUNTIME(hipcub::DeviceScan::ExclusiveSum(
    temp, temp_bytes,
    thrust::raw_pointer_cast(_neighbor_list->_d_max_neighbor_num.data()),
    thrust::raw_pointer_cast(_neighbor_list->_start_idx.data()),
    static_cast<int>(_linked_cell->_total_atoms_num)));
  CHECK_RUNTIME(FREE(temp));
  op::InitEndIndexOp<device::DEVICE_GPU> init_end_index_op;
  init_end_index_op(
      thrust::raw_pointer_cast(_neighbor_list->_d_neighbor_num.data()),
      thrust::raw_pointer_cast(_neighbor_list->_start_idx.data()),
      thrust::raw_pointer_cast(_neighbor_list->_end_idx.data()),
      _linked_cell->_total_atoms_num);
}