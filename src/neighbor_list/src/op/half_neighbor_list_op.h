#pragma once
#include "common/device_types.h"
#include "common/types.h"
#include "model/box.h"

namespace op {
template <typename DEVICE>
struct ComputeHalfNeighborsOp {
  void operator()(rbmd::Id* per_dimension_cells, rbmd::Id* neighbor_cell,
                  rbmd::Id neighbor_num, rbmd::Id total_cell,
                  rbmd::Id cell_count_within_cutoff);
};

// template <typename DEVICE>
// struct ComputeHalfNeighborsWithoutPBCOp {    //
// 这里生成邻居列表时要判断索引的
//   void operator()(rbmd::Id* neighbor_cell,
//                   rbmd::Id neighbor_num, rbmd::Id total_cell);
// };

template <typename DEVICE>
struct EstimateHalfNeighborListOp {
  void operator()(rbmd::Id* per_atom_cell_id,
                  rbmd::Id* in_atom_list_start_index,
                  rbmd::Id* in_atom_list_end_index, rbmd::Real cutoff_2,
                  rbmd::Id total_atom_num, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* neighbour_num,
                  rbmd::Id* max_neighbour_num, Box* box,
                  rbmd::Id* neighbor_cell, rbmd::Id neighbor_cell_num,
                  bool without_pbc);
};

template <typename DEVICE>
struct GenerateHalfNeighborListOp {
  void operator()(rbmd::Id* per_atom_cell_id,
                  rbmd::Id* in_atom_list_start_index,
                  rbmd::Id* in_atom_list_end_index, rbmd::Real cutoff_2,
                  rbmd::Id total_atom_num, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* max_neighbor_num,
                  rbmd::Id* neighbor_start, rbmd::Id* neighbor_end,
                  rbmd::Id* neighbors, Box* d_box, bool* should_realloc,
                  rbmd::Id* neighbor_cell, rbmd::Id neighbor_cell_num,
                  bool without_pbc);
};

template <>
struct ComputeHalfNeighborsOp<device::DEVICE_GPU> {
  void operator()(rbmd::Id* per_dimension_cells, rbmd::Id* neighbor_cell,
                  rbmd::Id neighbor_num, rbmd::Id total_cell,
                  rbmd::Id cell_count_within_cutoff);
};

template <>
struct EstimateHalfNeighborListOp<device::DEVICE_GPU> {
  void operator()(rbmd::Id* per_atom_cell_id,
                  rbmd::Id* in_atom_list_start_index,
                  rbmd::Id* in_atom_list_end_index, rbmd::Real cutoff_2,
                  rbmd::Id total_atom_num, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* neighbour_num,
                  rbmd::Id* max_neighbour_num, Box* box,
                  rbmd::Id* neighbor_cell, rbmd::Id neighbor_cell_num,
                  bool without_pbc);
};

template <>
struct GenerateHalfNeighborListOp<device::DEVICE_GPU> {
  void operator()(rbmd::Id* per_atom_cell_id,
                  rbmd::Id* in_atom_list_start_index,
                  rbmd::Id* in_atom_list_end_index, rbmd::Real cutoff_2,
                  rbmd::Id total_atom_num, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* max_neighbor_num,
                  rbmd::Id* neighbor_start, rbmd::Id* neighbor_end,
                  rbmd::Id* neighbors, Box* d_box, bool* should_realloc,
                  rbmd::Id* neighbor_cell, rbmd::Id neighbor_cell_num,
                  bool without_pbc);
};

}  // namespace op
