#include "common/device_types.h"
#include "common/rbmd_define.h"
#include "common/types.h"
#include "linked_cell/linked_cell.h"
#include "linked_cell_op.h"

namespace op {
__global__ void InitializeCell(LinkedCellDeviceDataPtr* linked_cell, Box* box,
                               Cell* cells, rbmd::Id total_cells) {
  __shared__ Int3 s_per_dimension_cells;
  __shared__ Real3 s_cell_length;
  __shared__ Real3 s_box_min;
  if (threadIdx.x == 0) {
    s_per_dimension_cells = make_Int3(linked_cell->_d_per_dimension_cells[0],
                                      linked_cell->_d_per_dimension_cells[1],
                                      linked_cell->_d_per_dimension_cells[2]);
    s_cell_length = make_Real3(linked_cell->_d_cell_length[0],
                               linked_cell->_d_cell_length[1],
                               linked_cell->_d_cell_length[2]);
    s_box_min =
        make_Real3(box->_coord_min[0], box->_coord_min[1], box->_coord_min[2]);
  }

  __syncthreads();

  unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < total_cells) {
    // One2threed
    int iz = idx / (s_per_dimension_cells.x * s_per_dimension_cells.y);
    int iy = (idx / s_per_dimension_cells.x) % s_per_dimension_cells.y;
    int ix = idx % s_per_dimension_cells.x;

    cells[idx]._cell_id = idx;

    cells[idx].cell_coord_min[0] = ix * s_cell_length.x + s_box_min.x;
    cells[idx].cell_coord_max[0] = (ix + 1) * s_cell_length.x + s_box_min.x;

    cells[idx].cell_coord_min[1] = iy * s_cell_length.y + s_box_min.y;
    cells[idx].cell_coord_max[1] = (iy + 1) * s_cell_length.y + s_box_min.y;

    cells[idx].cell_coord_min[2] = iz * s_cell_length.z + s_box_min.z;
    cells[idx].cell_coord_max[2] = (iz + 1) * s_cell_length.z + s_box_min.z;
  }
}

__global__ void AssignAtomsToCell(rbmd::Real* px, rbmd::Real* py,
                                  rbmd::Real* pz, Box* d_box,
                                  LinkedCellDeviceDataPtr* linked_cell,
                                  Cell* cells, rbmd::Id* per_atom_cell_id,
                                  rbmd::Id total_atoms_num) {
  unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= total_atoms_num) {
    return;
  }
  Int3 cell_idx;
  Real3 local_point = make_Real3(px[idx], py[idx], pz[idx]);
  // // TODO 可以取消这段逻辑吗？  这里local point 并没有修改原来的
  // if (local_point.x <= d_box->_coord_min[0]) {
  //   local_point.x += linked_cell->_d_cell_length[0] * 0.5;
  // } else if (local_point.x >= d_box->_coord_max[0]) {
  //   local_point.x -= linked_cell->_d_cell_length[0] * 0.5;
  // }
  //
  // if (local_point.y <= d_box->_coord_min[1]) {
  //   local_point.y += linked_cell->_d_cell_length[1] * 0.5;
  // } else if (local_point.y >= d_box->_coord_max[1]) {
  //   local_point.y -= linked_cell->_d_cell_length[1] * 0.5;
  // }
  //
  // if (local_point.z <= d_box->_coord_min[2]) {
  //   local_point.z += linked_cell->_d_cell_length[2] * 0.5;
  // } else if (local_point.z >= d_box->_coord_max[2]) {
  //   local_point.z -= linked_cell->_d_cell_length[2] * 0.5;
  // }

  cell_idx.x = MIN(MAX((int)(FLOOR((double)(local_point.x - d_box->_coord_min[0]) *
                                   linked_cell->_d_cell_length_reciprocal[0])) +
                           0,
                       0),
                   linked_cell->_d_per_dimension_cells[0] - 1);

  cell_idx.y = MIN(MAX((int)(FLOOR((double)(local_point.y - d_box->_coord_min[1]) *
                                   linked_cell->_d_cell_length_reciprocal[1])) +
                           0,
                       0),
                   linked_cell->_d_per_dimension_cells[1] - 1);

  cell_idx.z = MIN(MAX((int)(FLOOR((double)(local_point.z - d_box->_coord_min[2]) *
                                   linked_cell->_d_cell_length_reciprocal[2])) +
                           0,
                       0),
                   linked_cell->_d_per_dimension_cells[2] - 1);

  rbmd::Id cell_idx_1d =
      (cell_idx.z * linked_cell->_d_per_dimension_cells[1] + cell_idx.y) *
          linked_cell->_d_per_dimension_cells[0] +
      cell_idx.x;
  if (cells[cell_idx_1d].cell_coord_min[0] <= local_point.x &&
      cells[cell_idx_1d].cell_coord_min[1] <= local_point.y &&
      cells[cell_idx_1d].cell_coord_min[2] <= local_point.z &&
      local_point.x < cells[cell_idx_1d].cell_coord_max[0] &&
      local_point.y < cells[cell_idx_1d].cell_coord_max[1] &&
      local_point.z < cells[cell_idx_1d].cell_coord_max[2]) {
    per_atom_cell_id[idx] = cell_idx_1d;
    // 非特殊情况不会到下面这个分支
  }
  // else {
  //   if (local_point.x < cells[cell_idx_1d].cell_coord_min[0]) {
  //     cell_idx.x--;
  //   } else if (local_point.x >= cells[cell_idx_1d].cell_coord_max[0]) {
  //     cell_idx.x++;
  //   }
  //   if (local_point.y < cells[cell_idx_1d].cell_coord_min[1]) {
  //     cell_idx.y--;
  //   } else if (local_point.y >= cells[cell_idx_1d].cell_coord_max[1]) {
  //     cell_idx.y++;
  //   }
  //   if (local_point.z < cells[cell_idx_1d].cell_coord_min[2]) {
  //     cell_idx.z--;
  //   } else if (local_point.z >= cells[cell_idx_1d].cell_coord_max[2]) {
  //     cell_idx.z++;
  //   }
  //   cell_idx_1d =
  //       (cell_idx.z * linked_cell->_d_per_dimension_cells[1] + cell_idx.y) *
  //           linked_cell->_d_per_dimension_cells[0] +
  //       cell_idx.x;
  //   if (cells[cell_idx_1d].cell_coord_min[0] <= local_point.x &&
  //       cells[cell_idx_1d].cell_coord_min[1] <= local_point.y &&
  //       cells[cell_idx_1d].cell_coord_min[2] <= local_point.z &&
  //       local_point.x < cells[cell_idx_1d].cell_coord_max[0] &&
  //       local_point.y < cells[cell_idx_1d].cell_coord_max[1] &&
  //       local_point.z < cells[cell_idx_1d].cell_coord_max[2]) {
  //     per_atom_cell_id[idx] = cell_idx_1d;
  //   }
  // }
}

// 核函数来设置每个cell的start_index和end_index
/**
* @brief 计算每个单元格在排序后的原子列表中的起始和结束索引。
*
* 该函数通过并行化方式，扫描已排序的 `_per_atom_cell_id`（作为 `sorted_cell_index` 输入）
* 确定每个单元格对应的原子范围，并更新到 `d_in_atom_list_start_index` 和
* `d_in_atom_list_end_index` 中。
*
* 功能逻辑：
* - 对每个线程分配一个原子索引 `idx`，并根据其所属单元格ID进行以下操作：
*   1. 起始索引记录：当检测到当前原子是其单元格中第一个原子（首个原子或与前一原子的单元格ID不同），
*      更新 `d_in_atom_list_start_index[cell]` 为当前索引 `idx`。
*   2. 结束索引记录：当检测到当前原子是其单元格中最后一个原子（最后一个原子或与下一原子的单元格ID不同），
*      更新 `d_in_atom_list_end_index[cell]` 为当前索引的后一个值 `idx + 1`（区间为半开区间）。
*
* 输入参数：
* - `sorted_cell_index`：已按单元格ID排序的原子单元格ID数组，对应 `_per_atom_cell_id`。
* - `d_in_atom_list_start_index`：输出，每个单元格的起始索引，按单元格ID存储。
* - `d_in_atom_list_end_index`：输出，每个单元格的结束索引（半开区间）。
* - `num_atoms`：原子总数，表示数组的大小。
*
* 实现细节：
* - 使用 `atomicExch` 来保证在并行环境下对 `d_in_atom_list_start_index` 和
*   `d_in_atom_list_end_index` 的更新是线程安全的。性能上可能有优化空间。
* - `atomicExch` 操作：
*   - `d_in_atom_list_start_index[cell]` 在首次遇到单元格时更新。
*   - `d_in_atom_list_end_index[cell]` 在最后一次遇到单元格时更新。
*
*/
__global__ void ComputeCellRangesIndices(rbmd::Id* sorted_cell_index,
                                         rbmd::Id* d_in_atom_list_start_index,
                                         rbmd::Id* d_in_atom_list_end_index,
                                         rbmd::Id num_atoms) {
  unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx < num_atoms) {
    rbmd::Id cell = sorted_cell_index[idx];

    if (idx == 0 || cell != sorted_cell_index[idx - 1]) {
      //_d_in_atom_list_start_index[Cell] = idx;
      atomicExch(
          &d_in_atom_list_start_index[cell],
          idx);  // 一般来说不会同时设置，这个是后续性能优化的点，可以考虑不使用原子操作
    }

    if (idx == num_atoms - 1 || cell != sorted_cell_index[idx + 1]) {
      //_d_in_atom_list_end_index[Cell] = idx;
      atomicExch(&d_in_atom_list_end_index[cell],
                 idx + 1);  //  小于 不是小于等于
    }
  }
}

__global__ void MapAtomidToIdx(rbmd::Id* d_atomid2idx,
                               rbmd::Id* d_sorted_atom_idx,
                               rbmd::Id num_atoms) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx < num_atoms) {
    // 使用排序后的索引数组来创建映射
    d_atomid2idx[d_sorted_atom_idx[idx]] = idx;
  }
}

void InitializeCellOp<device::DEVICE_GPU>::operator()(
    LinkedCellDeviceDataPtr* linked_cell, Box* box, Cell* cells,
    rbmd::Id total_cells) {
  int threads_per_block = BLOCK_SIZE;
  int blocks_per_grid =
      (total_cells + threads_per_block - 1) / threads_per_block;
  CHECK_KERNEL(InitializeCell<<<blocks_per_grid, threads_per_block, 0, 0>>>(
      linked_cell, box, cells, total_cells));
}

void AssignAtomsToCellOp<device::DEVICE_GPU>::operator()(
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, Box* d_box,
    LinkedCellDeviceDataPtr* linked_cell, Cell* cells,
    rbmd::Id* per_atom_cell_id, rbmd::Id total_atoms_num) {
  unsigned int blocks_per_grid =
      (total_atoms_num + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(AssignAtomsToCell<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      px, py, pz, d_box, linked_cell, cells, per_atom_cell_id,
      total_atoms_num));
}

void ComputeCellRangesIndicesOp<device::DEVICE_GPU>::operator()(
    rbmd::Id* sorted_cell_index, rbmd::Id* d_in_atom_list_start_index,
    rbmd::Id* d_in_atom_list_end_index, rbmd::Id num_atoms) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(ComputeCellRangesIndices<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      sorted_cell_index, d_in_atom_list_start_index, d_in_atom_list_end_index,
      num_atoms));
}
void MapAtomidToIdxOp<device::DEVICE_GPU>::operator()(
    rbmd::Id* d_atomid2idx, rbmd::Id* d_sorted_atom_idx, rbmd::Id num_atoms) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(MapAtomidToIdx<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      d_atomid2idx, d_sorted_atom_idx, num_atoms));
}
}  // namespace op
