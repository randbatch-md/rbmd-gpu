#pragma once
#include <thrust/device_vector.h>

#include "common/types.h"

class NeighborList {
 public:
  explicit NeighborList(rbmd::Id total_atoms_num, bool is_half = false);

  bool* _h_is_half = nullptr;
  bool* _d_is_half = nullptr;
  thrust::device_vector<rbmd::Id> _d_neighbor_num{};
  thrust::device_vector<rbmd::Id> _d_max_neighbor_num{};
  rbmd::Id _h_total_max_neighbor_num = 0;
  thrust::device_vector<rbmd::Id> _d_neighbors{};
  thrust::device_vector<rbmd::Id> _start_idx{};
  thrust::device_vector<rbmd::Id> _end_idx{};
};
