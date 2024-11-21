#pragma once
#include "../neighbor_list/neighbor_list.h"
#include "base_neighbor_list_builder.h"

class FullNeighborListBuilder : public BaseNeighborListBuilder {
 public:
  explicit FullNeighborListBuilder();

  std::shared_ptr<NeighborList> Build() override;

 protected:
  /**
 * @brief 计算每个单元格的邻居单元格(包含自己)
 *
 * 此函数用于计算每个单元格的所有邻居单元格，并将其结果存储到 `_neighbor_cell` 向量中。
 * `_neighbor_cell` 的容量为 `total_cells * neighbor_cell_num`，用于保存每个单元格的邻居信息。
 *
 * 每个单元格的邻居可以通过以下方式访问：
 * `_neighbor_cell[cell_id * neighbor_cell_num + index]`，其中：
 * - `cell_id` 表示当前单元格的唯一 ID；
 * - `index` 的取值范围为 `[0, neighbor_cell_num)`。
 *
 * @note 此函数假定：
 * - 单元格数量 (`total_cells`) 和每个单元格的最大邻居数量 (`neighbor_cell_num`) 已预先定义；
 * - 周期性边界条件 (Periodic Boundary Conditions, PBC) 启用条件满足(总共的cell数量大于neighbor_cell_num)。
 *
 * 周期性边界条件的计算逻辑参考了以下文献：
 * Reich, S. *Numerical simulation in molecular dynamics: Numerics, algorithms, parallelization, applications*.
 * SIAM Review, 2010, 52(1): 213. （见第56页）
 */
  void ComputeNeighborCells() override;

  /**
  * @brief 计算每个单元格的邻居单元格(包含自己)
  * @note 不启用PBC,全部Cell均为邻居
  */
  void ComputeNeighborCellsWithoutPBC() override;

  /**
 * @brief 估计并调整邻居列表的大小。
 *
 * 该函数负责估计每个原子的最大邻居数量，并根据这些数量调整邻居列表数组的大小。
 * 具体步骤包括：
 * 1. 分配设备内存来存储总的最大的邻居数量。
 * 2. 将主机上的最大邻居数量复制到设备内存。
 * 3. 调用估计全邻居列表的操作，计算每个原子的邻居数量。
 * 4. 通过求和操作得到总的最大的邻居数量。
 * 5. 将设备内存中的总最大邻居数量复制回主机。
 * 6. 根据总的最大的邻居数量调整邻居列表数组的大小。
 * 7. 初始化邻居列表的索引。
 * 8. 设置 should_realloc 为 false，表示不需要重新分配内存。
 */
  void EstimateNeighborsList() override;

  rbmd::Id GenerateNeighborsList() override;

  std::shared_ptr<DeviceData> _device_data;
};
