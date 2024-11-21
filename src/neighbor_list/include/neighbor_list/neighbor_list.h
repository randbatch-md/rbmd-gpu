#pragma once
#include <thrust/device_vector.h>

#include "common/types.h"

/**
 * @class NeighborList
 * @brief 用于存储和管理原子邻居列表的数据结构，支持邻居列表的动态分配和操作。
 *
 * 该类非常适合Cuda及Rocm等平台，他完全一维存储了邻居列表信息。
 *
 * 该类为每个原子维护其邻居信息，包括邻居数量、最大容量、邻居索引等。
 * 提供了额外缓冲区以减少在模拟过程中因动态变化而导致的显存重新分配。
 */
class NeighborList {
 public:
  /**
   * @brief 构造函数
   * @param total_atoms_num 原子的总数，用于初始化数据结构。
   * @param is_half 是否启用半邻居表模式（默认关闭）。
   */
  explicit NeighborList(rbmd::Id total_atoms_num, bool is_half = false);
  /// @brief CPU 端标志指针，表示是否启用半邻居表模式。
  bool* _h_is_half = nullptr;
  /// @brief GPU 端标志指针，表示是否启用半邻居表模式。
  bool* _d_is_half = nullptr;
  /// 每个原子的邻居原子的实际数量（动态变化，表示当前邻居数）。
  thrust::device_vector<rbmd::Id> _d_neighbor_num{};
  /**
  * @brief 每个原子的最大邻居容量。
  *
  * 此值包括一定的缓冲空间，用于减少模拟过程中因邻居数变化导致的显存重新分配。
  */
  thrust::device_vector<rbmd::Id> _d_max_neighbor_num{};
  /// @brief 所有原子的最大邻居容量总和，用于分配邻居列表的总空间。
  rbmd::Id _h_total_max_neighbor_num = 0;
  /**
   * @brief 所有原子的邻居索引列表。
   *
   * 按以下方式排列：
   * - 原子0的邻居原子索引 + 缓冲区（最大容量）。
   * - 原子1的邻居原子索引 + 缓冲区（最大容量）。
   * - 以此类推，直到最后一个原子。
   *
   * 例如，若原子0的最大容量为 `N0`，实际邻居数为 `n0`，则其对应范围为 `[0, N0)`，其中 `[n0, N0)` 是缓冲区。
   */
  thrust::device_vector<rbmd::Id> _d_neighbors{};
  /// @brief 每个原子的邻居在 `_d_neighbors` 中的起始位置。
  thrust::device_vector<rbmd::Id> _start_idx{};
  /// @brief 每个原子的邻居在 `_d_neighbors` 中的结束位置（半开区间）。
  thrust::device_vector<rbmd::Id> _end_idx{};
  /**
  * @brief 每个原子的随机邻居索引列表，仅用于随机邻居选择（如 RBL 算法）。
  *
  * 数据排列方式为：
  * - 每个原子 `tid` 的随机邻居索引范围为 `[tid * neighbor_sample_num, (tid + 1) * neighbor_sample_num)`。
  * - 容量为 `原子数 * neighbor_sample_num`。
  */
  thrust::device_vector<rbmd::Id>
      _d_random_neighbor{};
  /**
  * @brief 每个原子随机挑选的邻居数量。
  *
  * 由于随机概率原因，实际挑选的邻居可能小于 `neighbor_sample_num`。
  */
  thrust::device_vector<rbmd::Id>
      _d_random_neighbor_num{};
  /// @brief 用于控制随机邻居选择频率的参数。
  rbmd::Id _selection_frequency = 0;
  // for i =_d_random_neighbor[tid*neighbor_sample_num]   i <
  // _d_random_neighbor_num[tid]

  // 输出邻居列表，便于调试
  void print(const std::string& filename);
};
