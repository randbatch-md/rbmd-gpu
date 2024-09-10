#pragma once
#include "common/types.h"
#include "model/box.h"

namespace op {
template <typename DEVICE>
struct ComputeFullNeighborsOp {
  void operator()(rbmd::Id* per_dimension_cells, rbmd::Id* neighbor_cell,
                  rbmd::Id neighbor_num, rbmd::Id total_cell,
                  rbmd::Id cell_count_within_cutoff);
};

template <typename DEVICE>
struct EstimateFullNeighborListOp {
  void operator()(rbmd::Id* per_atom_cell_id,
                  rbmd::Id* in_atom_list_start_index,
                  rbmd::Id* in_atom_list_end_index, rbmd::Real cutoff_2,
                  rbmd::Id total_atom_num, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* neighbour_num,
                  rbmd::Id* max_neighbour_num, Box* box,
                  rbmd::Id* neighbor_cell, rbmd::Id neighbor_cell_num);
};

template <typename DEVICE>
struct GenerateFullNeighborListOp {
  void operator()(rbmd::Id* per_atom_cell_id,
                  rbmd::Id* in_atom_list_start_index,
                  rbmd::Id* in_atom_list_end_index, rbmd::Real cutoff_2,
                  rbmd::Id total_atom_num, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* max_neighbor_num,
                  rbmd::Id* neighbor_start, rbmd::Id* neighbor_end,
                  rbmd::Id* neighbors, Box* d_box, bool* should_realloc,
                  rbmd::Id* neighbor_cell, rbmd::Id neighbor_cell_num);
};



template <>
struct ComputeFullNeighborsOp<device::DEVICE_GPU> {
  void operator()(rbmd::Id* per_dimension_cells, rbmd::Id* neighbor_cell,
                  rbmd::Id neighbor_num, rbmd::Id total_cell,
                  rbmd::Id cell_count_within_cutoff);
};

template <>
struct EstimateFullNeighborListOp<device::DEVICE_GPU> {
  void operator()(rbmd::Id* per_atom_cell_id,
                  rbmd::Id* in_atom_list_start_index,
                  rbmd::Id* in_atom_list_end_index, rbmd::Real cutoff_2,
                  rbmd::Id total_atom_num, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* neighbour_num,
                  rbmd::Id* max_neighbour_num, Box* box,
                  rbmd::Id* neighbor_cell, rbmd::Id neighbor_cell_num);
};

template <>
struct GenerateFullNeighborListOp<device::DEVICE_GPU> {
  void operator()(rbmd::Id* per_atom_cell_id,
                  rbmd::Id* in_atom_list_start_index,
                  rbmd::Id* in_atom_list_end_index, rbmd::Real cutoff_2,
                  rbmd::Id total_atom_num, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* max_neighbor_num,
                  rbmd::Id* neighbor_start, rbmd::Id* neighbor_end,
                  rbmd::Id* neighbors, Box* d_box, bool* should_realloc,
                  rbmd::Id* neighbor_cell, rbmd::Id neighbor_cell_num);
};

}  // namespace op
