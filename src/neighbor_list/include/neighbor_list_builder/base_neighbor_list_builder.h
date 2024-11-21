#pragma once
#include "../../data_manager/include/model/box.h"
#include "../common/object.h"
#include "../common/types.h"
#include "../neighbor_list/include/linked_cell/linked_cell.h"
#include "../neighbor_list/neighbor_list.h"

// linked Cell should new in out
class BaseNeighborListBuilder : public Object {
 public:
  explicit BaseNeighborListBuilder();
  ~BaseNeighborListBuilder() override;

  virtual std::shared_ptr<NeighborList> Build() = 0;

 protected:
  std::shared_ptr<LinkedCell> _linked_cell;
  std::shared_ptr<NeighborList> _neighbor_list = nullptr;

  virtual void ComputeNeighborCells() = 0;

  virtual void ComputeNeighborCellsWithoutPBC() = 0;

  virtual void EstimateNeighborsList() = 0;

  virtual rbmd::Id GenerateNeighborsList() = 0;

  // reduce求和
  void ReductionSum(rbmd::Id* d_src_array, rbmd::Id* d_dst, rbmd::Id size);

  /**
   * @brief 初始化邻居列表的开始和结束索引
   *
   * 此函数用于初始化每个原子的邻居列表开始和结束索引，并存储到相关数据结构中。
   * 邻居列表 (`neighborlist`) 是以一维数组的形式存储所有原子的邻居信息，
   * 因此需要通过扫描前缀和（Prefix Sum）计算开始索引 (`start_index`)，
   * 而结束索引由以下公式确定：
   *
   * end_index = start_index+ neighbor_num
   *
   * - `start_index`：表示当前原子在 `neighborlist` 中的起始位置；
   * - `neighbor_num`：表示当前原子的实际邻居数量；
   * - `end_index`：为邻居列表结束索引（不包括 `end_index` 本身）。
   *
   * @note 该函数假设每个原子的邻居数量 (`neighbor_num`) 已经预先计算。
   *
   * ### 示例
   * 假设有以下 3 个原子，其邻居数量分别为：
   * - 原子 0：邻居数量 = 2
   * - 原子 1：邻居数量 = 3
   * - 原子 2：邻居数量 = 1
   *
   * 邻居列表 `neighborlist` 的线性排列如下：
   * ```
   * neighborlist = [A, B, C, D, E, F]
   * ```
   * 计算结果：
   * - 原子 0：`start_index = 0`，`end_index = 0 + 2 = 2`，对应邻居为 `[A, B]`。
 * - 原子 1：`start_index = 2(0+2)`，`end_index = 2 + 3 = 5`，对应邻居为 `[C, D, E]`。
 * - 原子 2：`start_index = 5(0+2+3)`，`end_index = 5 + 1 = 6`，对应邻居为 `[F]`。
 */
  void InitNeighborListIndices();

  rbmd::Id _neighbor_cell_num = 0;
  rbmd::Id should_realloc = RBMD_TRUE;
  Box* _d_box;   // TODO 可能不太适合 待重构
  rbmd::Id* _d_should_realloc;
  rbmd::Real _trunc_distance_power_2 = 0;  //生成邻居的截断距离平方 通常为cutoff平方，rbl时为rcore平方

};
