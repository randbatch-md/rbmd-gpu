#pragma once
#include "common/device_types.h"
#include "common/types.h"
#include "model/box.h"

namespace op {
template <typename DEVICE>
struct ComputeMACENeighborsOp {
  void operator()(rbmd::Id* per_dimension_cells, rbmd::Id* neighbor_cell,
                  rbmd::Id neighbor_num, rbmd::Id total_cell,
                  rbmd::Id cell_count_within_cutoff);
};

template <typename DEVICE>
struct ComputeMACENeighborsWithoutPBCOp {
  void operator()(rbmd::Id* neighbor_cell, rbmd::Id neighbor_num,
                  rbmd::Id total_cell);
};

template <typename DEVICE>
struct EstimateMACENeighborListOp {
  void operator()(rbmd::Id* per_atom_cell_id,
                  rbmd::Id* in_atom_list_start_index,
                  rbmd::Id* in_atom_list_end_index, rbmd::Real cutoff_2,
                  rbmd::Id total_atom_num, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* neighbour_num,
                  rbmd::Id* max_neighbour_num, Box box,
                  rbmd::Id* neighbor_cell, rbmd::Id neighbor_cell_num);
};

template <typename DEVICE>
struct GenerateMACENeighborListOp {
  void operator()(rbmd::Id* per_atom_cell_id,
                  rbmd::Id* in_atom_list_start_index,
                  rbmd::Id* in_atom_list_end_index, rbmd::Real cutoff_2,
                  rbmd::Id total_atom_num, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* max_neighbor_num,
                  rbmd::Id* neighbor_start, rbmd::Id* neighbor_end,
                  rbmd::Id* neighbors, Box box, rbmd::Id* should_realloc,
                  rbmd::Id* neighbor_cell, rbmd::Id neighbor_cell_num);
};

template <>
struct ComputeMACENeighborsOp<device::DEVICE_GPU> {
  void operator()(rbmd::Id* per_dimension_cells, rbmd::Id* neighbor_cell,
                  rbmd::Id neighbor_num, rbmd::Id total_cell,
                  rbmd::Id cell_count_within_cutoff);
};

template <>
struct ComputeMACENeighborsWithoutPBCOp<device::DEVICE_GPU> {
  void operator()(rbmd::Id* neighbor_cell, rbmd::Id neighbor_num,
                  rbmd::Id total_cell);
};

template <>
struct EstimateMACENeighborListOp<device::DEVICE_GPU> {
  void operator()(rbmd::Id* per_atom_cell_id,
                  rbmd::Id* in_atom_list_start_index,
                  rbmd::Id* in_atom_list_end_index, rbmd::Real cutoff_2,
                  rbmd::Id total_atom_num, rbmd::Real* px, rbmd::Real* py,
                  rbmd::Real* pz, rbmd::Id* neighbour_num,
                  rbmd::Id* max_neighbour_num, Box box,
                  rbmd::Id* neighbor_cell, rbmd::Id neighbor_cell_num);
};

template <>
struct GenerateMACENeighborListOp<device::DEVICE_GPU> {
  void operator()(    rbmd::Id*  per_atom_cell_id,
    rbmd::Id*  in_atom_list_start_index,
    rbmd::Id*  in_atom_list_end_index, rbmd::Real cutoff_2,
    rbmd::Id total_atom_num,  rbmd::Real*  px,
    rbmd::Real*  py, rbmd::Real*  pz,
    rbmd::Id*  max_neighbor_num,
    rbmd::Id*  neighbor_start, rbmd::Id*  neighbor_end,
    rbmd::Id*  neighbors, Box  d_box,
        rbmd::Id*  neighbors_atoms, rbmd::Real*   unit_shiftx,
        rbmd::Real*   unit_shifty,  rbmd::Real*  unit_shiftz,
        rbmd::Real*   shiftx,	rbmd::Real*   shifty,
        rbmd::Real*   shiftz,
    rbmd::Id*  should_realloc, rbmd::Id*  neighbor_cell,
    rbmd::Id neighbor_cell_num);
};

}  // namespace op

__host__ __device__ __forceinline__ rbmd::Real CaculateDistance(
    Box box, rbmd::Real i_x, rbmd::Real i_y, rbmd::Real i_z, rbmd::Real j_x,
    rbmd::Real j_y, rbmd::Real j_z) {
  rbmd::Real dx = i_x - j_x;
  rbmd::Real dy = i_y - j_y;
  rbmd::Real dz = i_z - j_z;
  MinImageDistance(box, dx, dy, dz);
  return (POW(dx, 2) + POW(dy, 2) + POW(dz, 2));
}