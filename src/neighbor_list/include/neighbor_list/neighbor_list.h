#pragma once
#include <thrust/device_vector.h>

#include "common/types.h"

class NeighborList {
 public:
  explicit NeighborList(rbmd::Id total_atoms_num, bool is_half = false);

  bool* _h_is_half = nullptr;
  bool* _d_is_half = nullptr;
  thrust::device_vector<rbmd::Id> _d_neighbor_num{};  // 每个原子的邻居数量
  thrust::device_vector<rbmd::Id> _d_max_neighbor_num{};  // 每个原子的最大邻居容量
  rbmd::Id _h_total_max_neighbor_num = 0;
  thrust::device_vector<rbmd::Id> _d_neighbors{}; // 每个原子的邻居列表  when use rbl is: s_neighbor
  thrust::device_vector<rbmd::Id> _start_idx{}; // 每个原子的邻居在_d_neighbors的起始位置
  thrust::device_vector<rbmd::Id> _end_idx{}; // 每个原子的邻居在_d_neighbors的结束位置
  thrust::device_vector<rbmd::Id> _d_random_neighbor{}; // just for rbl  cs neighbor  开始：atom_id  * neighbor_sample_num
  thrust::device_vector<rbmd::Id> _d_random_neighbor_num{}; // 概率原因，每个原子的随机邻居可能不到neighbor_sample_num个，这里记录每个原子随机挑选的邻居个数
  rbmd::Id _selection_frequency = 0;
  // for i =_d_random_neighbor[tid*neighbor_sample_num]   i <  _d_random_neighbor_num[tid]
  void print(const std::string& filename);
};
