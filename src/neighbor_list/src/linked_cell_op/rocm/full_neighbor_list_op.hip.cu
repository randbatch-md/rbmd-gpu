#include <hip/hip_runtime.h>

#include <hipcub/hipcub.hpp>

#include "common/device_types.h"
#include "common/rbmd_define.h"
#include "common/types.h"
#include "full_neighbor_list_op.h"
#include "model/box.h"

namespace op {
__global__ void ComputeFullNeighbors(rbmd::Id* per_dimension_cells,
                                     rbmd::Id* neighbor_cell,
                                     rbmd::Id neighbor_num, rbmd::Id total_cell,
                                     rbmd::Id cell_count_within_cutoff) {
  const unsigned int cell_idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (cell_idx < total_cell) {
    rbmd::Id idx_z =
        cell_idx / (per_dimension_cells[0] * per_dimension_cells[1]);
    rbmd::Id idx_y =
        (cell_idx - idx_z * per_dimension_cells[0] * per_dimension_cells[1]) /
        per_dimension_cells[0];
    rbmd::Id idx_x = cell_idx - per_dimension_cells[0] *
                                    (idx_y + per_dimension_cells[1] * idx_z);
    int neighbor_count = 0;
    for (int dz = -cell_count_within_cutoff; dz <= cell_count_within_cutoff;
         ++dz) {
      // TODO cell in cutoff
      // RBL 2！
      for (int dy = -cell_count_within_cutoff; dy <= cell_count_within_cutoff;
           ++dy) {
        for (int dx = -cell_count_within_cutoff; dx <= cell_count_within_cutoff;
             ++dx) {
          int neighbor_z =
              (idx_z + dz + per_dimension_cells[2]) % per_dimension_cells[2];
          int neighbor_y =
              (idx_y + dy + per_dimension_cells[1]) % per_dimension_cells[1];
          int neighbor_x =
              (idx_x + dx + per_dimension_cells[0]) % per_dimension_cells[0];

          rbmd::Id neighbour_cell_idx =
              (neighbor_z * per_dimension_cells[1] + neighbor_y) *
                  per_dimension_cells[0] +
              neighbor_x;

          neighbor_cell[cell_idx * neighbor_num + neighbor_count] =
              neighbour_cell_idx;
          ++neighbor_count;
        }
      }
    }
  }
}

// use while total_cell < computed neighbor_num
__global__ void ComputeFullNeighborsWithoutPBC(rbmd::Id* neighbor_cell, rbmd::Id neighbor_num,
                                         rbmd::Id total_cell) {
  const unsigned int cell_idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (cell_idx < total_cell) {
    int neighbor_count = 0;
    for (int neighbor_cell_idx = 0; neighbor_cell_idx < neighbor_num; ++neighbor_cell_idx) {
      neighbor_cell[cell_idx * neighbor_num + neighbor_count] =
        neighbor_cell_idx;
      ++neighbor_count;
    }
  }
}


__host__ __device__ __forceinline__ rbmd::Real CaculateDistance(
    Box* box, rbmd::Real i_x, rbmd::Real i_y, rbmd::Real i_z, rbmd::Real j_x,
    rbmd::Real j_y, rbmd::Real j_z) {
  rbmd::Real dx = i_x - j_x;
  rbmd::Real dy = i_y - j_y;
  rbmd::Real dz = i_z - j_z;
  MinImageDistance(box, dx, dy, dz);
  return (POW(dx, 2) + POW(dy, 2) + POW(dz, 2));
}

/// 一个原子对应一个线程束   建议调用的时候 -EPSILON
__global__ void EstimateFullNeighborList(
    rbmd::Id* per_atom_cell_id, rbmd::Id* in_atom_list_start_index,
    rbmd::Id* in_atom_list_end_index, rbmd::Real cutoff_2,
    rbmd::Id total_atom_num, rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz,
    rbmd::Id* neighbour_num, rbmd::Id* max_neighbour_num, Box* box,
    rbmd::Id* neighbor_cell, rbmd::Id neighbor_cell_num) {
  // cutoff2是平方
  extern __shared__ hipcub::WarpReduce<rbmd::Id>::TempStorage
      reduce_temp_storage[];
  rbmd::Id atom_neighbor_num = MIN_NBNUM;

  // 计算当前线程处理的原子的索引    ---- 当前线程原子就是atom_idx
  unsigned int atom_idx = (blockIdx.x * blockDim.x + threadIdx.x) / warpSize;
  if (atom_idx < total_atom_num) {
    rbmd::Real distance = 0;
    // 计算当前线程在warp中的位置
    unsigned int lane_id = (blockIdx.x * blockDim.x + threadIdx.x) % warpSize;
    atom_neighbor_num = 0;
    // 当前所属的cell
    rbmd::Id cell_idx = per_atom_cell_id[atom_idx];
    for (int i = 0; i < neighbor_cell_num; ++i) {
      rbmd::Id neighbour_cell_idx =
          neighbor_cell[cell_idx * neighbor_cell_num + i];
      for (rbmd::Id neighbor_atom_idx =
               in_atom_list_start_index[neighbour_cell_idx] + lane_id;
           neighbor_atom_idx < in_atom_list_end_index[neighbour_cell_idx];
           neighbor_atom_idx += warpSize) {
        if (atom_idx != neighbor_atom_idx) {
          distance =
              CaculateDistance(box, px[atom_idx], py[atom_idx], pz[atom_idx],
                               px[neighbor_atom_idx], py[neighbor_atom_idx],
                               /// TODO
                               /// 后续考虑是否要将数据先进一遍寄存器
                               pz[neighbor_atom_idx]);
          if (distance < cutoff_2) {
            atom_neighbor_num++;
          }
        }
      }
    }
    atom_neighbor_num = hipcub::WarpReduce<rbmd::Id>(
                            reduce_temp_storage[threadIdx.x / warpSize])
                            .Sum(atom_neighbor_num);
    if (lane_id == 0) {
      neighbour_num[atom_idx] = atom_neighbor_num;
      // 线程束对齐
      max_neighbour_num[atom_idx] =
          MAX(((rbmd::Id)ceil(atom_neighbor_num * SAFE_ZONE) + warpSize - 1) /
                  warpSize * warpSize,
              MIN_NBNUM);
    }
  }
}

__global__ void GenerateFullNeighborList(
    rbmd::Id* per_atom_cell_id, rbmd::Id* in_atom_list_start_index,
    rbmd::Id* in_atom_list_end_index, rbmd::Real cutoff_2,
    rbmd::Id total_atom_num, rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz,
    rbmd::Id* max_neighbor_num, rbmd::Id* neighbor_start,
    rbmd::Id* neighbor_end, rbmd::Id* neighbors, Box* d_box,
    bool* should_realloc, rbmd::Id* neighbor_cell, rbmd::Id neighbor_cell_num) {
  const unsigned int atom_idx =
      (blockIdx.x * blockDim.x + threadIdx.x) / warpSize;
  const unsigned int lane_id =
      (blockIdx.x * blockDim.x + threadIdx.x) % warpSize;

  __shared__ hipcub::WarpScan<int>::TempStorage temp_storage[WARP_SIZE];

  rbmd::Id neighbor_num = 0;

  if (atom_idx < total_atom_num) {
    neighbor_num = neighbor_start[atom_idx];
    rbmd::Id cell_idx = per_atom_cell_id[atom_idx];

    for (int i = 0; i < neighbor_cell_num; ++i) {
      rbmd::Id neighbour_cell_idx =
          neighbor_cell[cell_idx * neighbor_cell_num + i];
      rbmd::Id start = in_atom_list_start_index[neighbour_cell_idx];
      rbmd::Id end = in_atom_list_end_index[neighbour_cell_idx];

      for (rbmd::Id base_idx = start; base_idx < end; base_idx += warpSize) {
        rbmd::Id neighbor_atom_idx = base_idx + lane_id;
        int is_neighbor = 0;

        if (neighbor_atom_idx < end && atom_idx != neighbor_atom_idx) {
          rbmd::Real distance =
              CaculateDistance(d_box, px[atom_idx], py[atom_idx], pz[atom_idx],
                               px[neighbor_atom_idx], py[neighbor_atom_idx],
                               pz[neighbor_atom_idx]);

          if (distance < cutoff_2) {
            is_neighbor = 1;
          }
        }

        int offset;
        hipcub::WarpScan<int>(temp_storage[threadIdx.x / warpSize])
            .ExclusiveSum(is_neighbor, offset);
        if (is_neighbor) {
          neighbors[neighbor_num + offset] = neighbor_atom_idx;
        }

        neighbor_num += hipcub::ShuffleIndex<WARP_SIZE, rbmd::Id>(
            offset + is_neighbor, warpSize - 1,
            0xffffffff);  // 广播线程束最后一个
        // __shfl_sync(0xffffffff, offset + is_neighbor, warpSize - 1);
        // 兼容dcu
      }
    }

    if (lane_id == 0) {
      neighbor_end[atom_idx] = neighbor_num;
      rbmd::Id my_total_neighbor_num = neighbor_num - neighbor_start[atom_idx];
      if (my_total_neighbor_num > max_neighbor_num[atom_idx]) {
        *should_realloc = true;
      }
    }
  }
}

void ComputeFullNeighborsOp<device::DEVICE_GPU>::operator()(
    rbmd::Id* per_dimension_cells, rbmd::Id* neighbor_cell,
    rbmd::Id neighbor_num, rbmd::Id total_cell,
    rbmd::Id cell_count_within_cutoff) {
  unsigned int blocks_per_grid = (total_cell + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(ComputeFullNeighbors<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      per_dimension_cells, neighbor_cell, neighbor_num, total_cell,
      cell_count_within_cutoff));
}

void ComputeFullNeighborsWithoutPBCOp<device::DEVICE_GPU>::operator()(rbmd::Id* neighbor_cell, rbmd::Id neighbor_num, rbmd::Id total_cell){
  unsigned int blocks_per_grid = (total_cell + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(ComputeFullNeighborsWithoutPBC<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(neighbor_cell, neighbor_num, total_cell));
}

void EstimateFullNeighborListOp<device::DEVICE_GPU>::operator()(
    rbmd::Id* per_atom_cell_id, rbmd::Id* in_atom_list_start_index,
    rbmd::Id* in_atom_list_end_index, rbmd::Real cutoff_2,
    rbmd::Id total_atom_num, rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz,
    rbmd::Id* neighbour_num, rbmd::Id* max_neighbour_num, Box* box,
    rbmd::Id* neighbor_cell, rbmd::Id neighbor_cell_num) {
  unsigned int threads_per_block = BLOCK_SIZE;
  unsigned int warps_per_block = threads_per_block / WARP_SIZE;
  unsigned int blocks_per_grid =
      (total_atom_num + warps_per_block - 1) / warps_per_block;
  CHECK_KERNEL(
      EstimateFullNeighborList<<<blocks_per_grid, BLOCK_SIZE,
                                 sizeof(hipcub::WarpReduce<int>::TempStorage) *
                                     (BLOCK_SIZE / WARP_SIZE),
                                 0>>>(
          per_atom_cell_id, in_atom_list_start_index, in_atom_list_end_index,
          cutoff_2, total_atom_num, px, py, pz, neighbour_num,
          max_neighbour_num, box, neighbor_cell, neighbor_cell_num));
}

void GenerateFullNeighborListOp<device::DEVICE_GPU>::operator()(
    rbmd::Id* per_atom_cell_id, rbmd::Id* in_atom_list_start_index,
    rbmd::Id* in_atom_list_end_index, rbmd::Real cutoff_2,
    rbmd::Id total_atom_num, rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz,
    rbmd::Id* max_neighbor_num, rbmd::Id* neighbor_start,
    rbmd::Id* neighbor_end, rbmd::Id* neighbors, Box* d_box,
    bool* should_realloc, rbmd::Id* neighbor_cell, rbmd::Id neighbor_cell_num) {
  unsigned int threads_per_block = BLOCK_SIZE;
  unsigned int warps_per_block = threads_per_block / WARP_SIZE;
  unsigned int blocks_per_grid =
      (total_atom_num + warps_per_block - 1) / warps_per_block;
  CHECK_KERNEL(
      GenerateFullNeighborList<<<blocks_per_grid, BLOCK_SIZE,
                                 sizeof(hipcub::WarpScan<int>::TempStorage) *
                                     (BLOCK_SIZE / WARP_SIZE),
                                 0>>>(
          per_atom_cell_id, in_atom_list_start_index, in_atom_list_end_index,
          cutoff_2, total_atom_num, px, py, pz, max_neighbor_num,
          neighbor_start, neighbor_end, neighbors, d_box, should_realloc,
          neighbor_cell, neighbor_cell_num));

  }

} // namespace op
