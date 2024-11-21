#pragma once
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include "../../data_manager/include/model/box.h"
#include "cell.h"
#include "common/rbmd_define.h"
#include "common/types.h"
#include "config_data.h"
#include "model/device_data.h"
#include "model/structure_info_data.h"

/**
 * @brief 用于方便地在核函数中获取所需使用的LinkedCell相关的参数
 * 参数:
 * -_d_total_cells 总共的cell数量
 * -_d_per_dimension_cells 盒子范围内每个维度划分的cell数量
 * -_d_cell_length 每个cell的x,y,z长度(浮点数表示)
 * -_d_cutoff cutoff参数，可参考：Rapaport D C. The art of molecular dynamics
 * simulation[M]. Cambridge university press, 2004.
 * -_d_cell_length_reciprocal cell长度的倒数，方便计算的同时节省了性能
 *
 * ALIGN：用于内存对齐
 *
 */
struct LinkedCellDeviceDataPtr {
  rbmd::Id _d_total_cells = 0;
  rbmd::Id ALIGN(ALIGN_SIZE(rbmd::Id, 3)) _d_per_dimension_cells[3]{};
  rbmd::Real ALIGN(ALIGN_SIZE(rbmd::Real, 3)) _d_cell_length[3]{};
  rbmd::Real _d_cutoff = 0;
  // Cell* _d_cells = nullptr;
  rbmd::Real ALIGN(ALIGN_SIZE(rbmd::Real, 3)) _d_cell_length_reciprocal[3]{};
};

/**
 * @brief LinkedCell类用于实现分子动力学模拟中的链表单元格方法（Linked Cell
 * Method）。 该方法是一种空间划分技术，用于优化分子间相互作用的计算效率。
 *
 * 链表单元格方法的基本思想是将模拟空间划分为多个小的单元格，并将每个粒子分配到相应的单元格中。
 * 然后，只计算位于同一单元格或相邻单元格中的粒子之间的相互作用，从而减少不必要的计算。
 *
 * 参考文献：
 * - Reich S. Numerical simulation in molecular dynamics: Numerics, algorithms,
 * parallelization, applications[J]. SIAM Review, 2010, 52(1): 213.
 *
 * 主要功能：
 * 1. 将模拟空间划分为单元格。
 * 2. 将粒子分配到相应的单元格中。
 * 3. 根据原子所属的cell重新排序
 *
 */
class LinkedCell {
 public:
  LinkedCell();
  ~LinkedCell();

  /// 总共的cell数量
  rbmd::Id _total_cells = 0;
  /// 总共的原子数量
  rbmd::Id _total_atoms_num = 0;
  /// 各维度的cell的数量
  rbmd::Id ALIGN(ALIGN_SIZE(rbmd::Id, 3)) _per_dimension_cells[3]{};
  /// 各维度cell的长度
  rbmd::Real ALIGN(ALIGN_SIZE(rbmd::Real, 3)) _cell_length[3]{};
  /// 各维度cell的长度的倒数
  rbmd::Real ALIGN(ALIGN_SIZE(rbmd::Real, 3)) _cell_length_reciprocal[3]{};

  /// 每个cell内的原子在原子列表的开始索引
  thrust::device_vector<rbmd::Id> _in_atom_list_start_index{};
  /// 每个cell内的原子在原子列表的结束索引
  thrust::device_vector<rbmd::Id> _in_atom_list_end_index{};
  /// 所有cell 行主序（Row-Major Order) 存储。
  /// 在内存中，先存储 x 方向上的单元格，然后是 y 方向上的单元格，最后是 z
  /// 方向上的单元格。也就是说，x 方向上的变化最快，y 方向上的变化次之，z
  /// 方向上的变化最慢
  thrust::device_vector<Cell> _cells{};
  /// 每个cell的邻居cell(通常cell_length =
  /// cutoff,此情况下如果总的cell超过27个则每个cell有27个邻居)
  thrust::device_vector<rbmd::Id> _neighbor_cell{};
  /// 由于SortAtomsByCellKey对原子列表进行了排序，在计算键角等需要使用id来进行查找位置时需要映射
  /// TODO 统一处理键角的Table
  thrust::device_vector<rbmd::Id> _atom_id_to_idx{};
  /// 模拟参数 非键相互作用的截断距离（cutoff
  /// distance）。这个参数决定了在计算原子间的非键相互作用（如范德华力和静电相互作用）时，哪些原子对会被考虑。单位通常为埃
  /// (Å)
  rbmd::Real _cutoff = 0;
  /// 一个cutoff长度内划分cell的个数。默认为1，即cutoff=cell_length
  rbmd::Id _cell_count_within_cutoff = 1;
  /// 每个原子所属的cell的id
  thrust::device_vector<rbmd::Id> _per_atom_cell_id{};

  /**
   * @brief 构建给定模拟盒的链表单元结构。
   *
   * 该函数根据截断距离和截断范围内的单元数量将模拟盒划分为单元。它计算每个维度中的单元数量、
   * 总单元数量以及每个单元的长度。它还初始化用于存储每个单元内原子索引的必要数据结构。
   *
   * @param[in,out] box 要划分为单元的模拟盒。
   *
   * @note 该函数旨在主机（CPU）端调用。
   */
  __host__ void Build(Box* box);

  /**
   * @brief 获取指向 LinkedCellDeviceDataPtr 结构的指针。
   *
   * 该函数返回一个指向 LinkedCellDeviceDataPtr
   * 结构的指针，该结构包含了在核函数中使用的所有 LinkedCell 相关参数。
   * 这些参数包括总单元数量、每个维度划分的单元数量、每个单元的长度、截断参数以及单元长度的倒数。
   *
   * @return LinkedCellDeviceDataPtr* 指向 LinkedCellDeviceDataPtr 结构的指针。
   *
   * @note 该函数旨在主机（CPU）端调用，返回的指针指向设备（GPU）上的数据结构。
   */
  LinkedCellDeviceDataPtr* GetDataPtr();

  /**
   * @brief 初始化单元格，为每个单元格设置ID及左上角和右下角的坐标。
   *
   * 该函数遍历所有的单元格，并为每个单元格设置唯一的ID，同时计算并设置
   * 单元格的左上角和右下角的坐标。
   */
  void InitializeCells();

  /**
   * @brief 根据原子坐标为每个原子分配所属的单元格ID，并存储在
   * `_per_atom_cell_id` 中。
   *
   * 该函数遍历所有原子，根据其坐标计算出所属的单元格，并将该单元格的ID存储在
   * `_per_atom_cell_id` 中。用于存储每个原子所属的单元格ID。
 */
  void AssignAtomsToCell();

  /**
   * @brief 对基于单元格键值的原子列表进行排序，并更新原子ID到索引的映射关系。
   *
   * 此函数用于将多个与原子相关的一维数组（SoA形式）组合，并按 `_per_atom_cell_id`
   * 中的单元格ID对这些数组进行稳定排序。排序后的结果是：
   * - 原子的排列顺序根据它们所属的单元格ID升序。
   * - `_per_atom_cell_id` 中的键值（单元格ID）也按升序排列（参考 `stable_sort_by_key` 的原理）。
   *
   * 排序结果：
   * - 排列后的SoA：`cell0` 的所有原子 - `cell1` 的所有原子 - ...（按单元格顺序）。
   * - 排列后的 `_per_atom_cell_id`：`cell0` 的所有键值（`0*n0`） - `cell1` 的所有键值（`1*n1`） - ...
   *
   * 排序完成后，函数还会执行 `MapAtomidToIdxOp` 操作，**更新原子ID到索引（Idx）的映射**，
   * 以便后续模拟中，除键角处理外，所有操作均基于索引进行。
   *
   * 参考文献：
   * Kylasa, Sudhir B.; Aktulga, Hasan; and Grama, Ananth, "PG-PuReMD: A Parallel-GPU Reactive Molecular
   * Dynamics Package" (2013). Department of Computer Science Technical Reports. Paper 1768.
   * https://docs.lib.purdue.edu/cstech/1768
   *
   */
  void SortAtomsByCellKey();

  /**
   * @brief 根据 `_per_atom_cell_id` 计算每个单元格在原子列表中的范围。
   *
   * 由于原子列表已经按照单元格ID进行排序，该函数可以很容易地计算出每个单元格
   * 包含多少个原子，以及它在原子列表中的开始和结束索引。这些信息将用于后续的
   * 邻近搜索和处理。
   *
   * 参考：Kylasa, Sudhir B.; Aktulga, Hasan; and Grama, Ananth, "PG-PuReMD: A Parallel-GPU Reactive Molecular
   * Dynamics Package" (2013).Department of Computer Science Technical Reports. Paper 1768.
   * https://docs.lib.purdue.edu/cstech/1768
   */
  void ComputeCellRangesIndices();

  /**
   * @brief 内存拷贝的操作，会自动触发。
   */
  void SyncHToD();


  /**
   * 由于SortAtomsByCellKey做了映射，读取文件计算键角力(在邻居力计算之后)需要将Id进行一次映射。调用这个函数。
   * @param d_target 需要映射的原子ID的数组
   * TODO 并行考虑范围问题？
   */
  template <typename T>
  void MapAtomId(thrust::device_vector<T>& d_target);

  /**
   * 由于SortAtomsByCellKey做了映射,如果计算键力角力在构建邻居列表之前(或同时)，这力加和时需要调用这个函数
   * @param d_target 键角力数组
   */
  template <typename T>
 void MapAtomIdForce(thrust::device_vector<T>& d_target);

 private:
  LinkedCellDeviceDataPtr* _linked_cell_device_data_ptr = nullptr;
  std::shared_ptr<DeviceData> _device_data;
  std::shared_ptr<StructureInfoData> _structure_info_data;
  std::shared_ptr<ConfigData> _config_data;
  /**
   * @brief 分配设备显存，会自动触发。
   */
  void AllocDeviceMemory();
};
