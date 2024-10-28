#pragma once
#include "common/device_types.h"
#include "common/types.h"
#include "model/box.h"
namespace op {
template <typename DEVICE>
struct GenerateRblFullNeighborListOp {
  void operator()(rbmd::Id* per_atom_cell_id,
                  rbmd::Id* in_atom_list_start_index,
                  rbmd::Id* in_atom_list_end_index,
                  rbmd::Real trunc_distance_power_2, rbmd::Real cutoff_2,
                  rbmd::Id total_atom_num, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* max_neighbor_num,
                  rbmd::Id* neighbor_start, rbmd::Id* neighbor_end,
                  rbmd::Id* neighbors, rbmd::Id neighbor_sample_num,
                  rbmd::Id* random_neighbors, rbmd::Id* random_neighbors_num,
                  Box* d_box, rbmd::Id* should_realloc, rbmd::Id* neighbor_cell,
                  rbmd::Id neighbor_cell_num,rbmd::Id selection_frequency);
};

template <>
struct GenerateRblFullNeighborListOp<device::DEVICE_GPU> {
  void operator()(rbmd::Id* per_atom_cell_id,
                  rbmd::Id* in_atom_list_start_index,
                  rbmd::Id* in_atom_list_end_index,
                  rbmd::Real trunc_distance_power_2, rbmd::Real cutoff_2,
                  rbmd::Id total_atom_num, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* max_neighbor_num,
                  rbmd::Id* neighbor_start, rbmd::Id* neighbor_end,
                  rbmd::Id* neighbors, rbmd::Id neighbor_sample_num,
                  rbmd::Id* random_neighbors, rbmd::Id* random_neighbors_num,
                  Box* d_box, rbmd::Id* should_realloc, rbmd::Id* neighbor_cell,
                  rbmd::Id neighbor_cell_num,rbmd::Id selection_frequency);
};

}  // namespace op