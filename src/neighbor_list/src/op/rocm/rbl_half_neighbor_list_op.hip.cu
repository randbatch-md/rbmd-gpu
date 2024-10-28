#include "full_neighbor_list_op.h"
#include "rbl_half_neighbor_list_op.h"
#include "src/op/rbl_half_neighbor_list_op.h"
namespace op {

// TODO warp op and optimize
//! Note 暂时未测试Warp的
//! 因为是用全部邻居算，所以和半邻居的不一样  大部分代码重复 需要重构
__global__ void GenerateRBLHalfNeighborList(
    rbmd::Id* __restrict__  per_atom_cell_id, rbmd::Id* __restrict__  in_atom_list_start_index,
    rbmd::Id* __restrict__  in_atom_list_end_index, rbmd::Real trunc_distance_power_2,
    rbmd::Real cutoff_2, rbmd::Id total_atom_num, rbmd::Real* __restrict__  px,
    rbmd::Real* __restrict__  py, rbmd::Real* __restrict__  pz, rbmd::Id* __restrict__  max_neighbor_num,
    rbmd::Id* __restrict__  neighbor_start, rbmd::Id* __restrict__  neighbor_end, rbmd::Id* __restrict__  neighbors,
    rbmd::Id neighbor_sample_num, rbmd::Id* __restrict__  random_neighbors,
    rbmd::Id* __restrict__  random_neighbors_num, Box* __restrict__  d_box, rbmd::Id* __restrict__  should_realloc,
    rbmd::Id* __restrict__  neighbor_cell, rbmd::Id neighbor_cell_num,
    rbmd::Id selection_frequency) {
  const unsigned int atom_idx = (blockIdx.x * blockDim.x + threadIdx.x);
  rbmd::Real distance = 0;
  rbmd::Id neighbor_num = 0;
  rbmd::Id random_neighbor_num = 0;
  rbmd::Id index_shell = 0;
  rbmd::Id random_neighbor_start = 0;
  __shared__ rbmd::Real shared_px[BLOCK_SIZE];
  __shared__ rbmd::Real shared_py[BLOCK_SIZE];
  __shared__ rbmd::Real shared_pz[BLOCK_SIZE];
  if (atom_idx < total_atom_num) {
    neighbor_num = neighbor_start[atom_idx];
    random_neighbor_start = atom_idx * neighbor_sample_num;
    rbmd::Id cell_idx = __ldg(&per_atom_cell_id[atom_idx]);

    // 缓存原子位置到共享内存
    if (threadIdx.x < BLOCK_SIZE) {
      shared_px[threadIdx.x] = __ldg(&px[atom_idx]);
      shared_py[threadIdx.x] = __ldg(&py[atom_idx]);
      shared_pz[threadIdx.x] = __ldg(&pz[atom_idx]);
    }
    __syncthreads();
    for (int i = 0; i < neighbor_cell_num; ++i) {
      rbmd::Id neighbour_cell_idx =
               __ldg(&neighbor_cell[cell_idx * neighbor_cell_num + i]);
      rbmd::Id start = __ldg(&in_atom_list_start_index[neighbour_cell_idx]);
      rbmd::Id end = __ldg(&in_atom_list_end_index[neighbour_cell_idx]);
      for (rbmd::Id neighbor_atom_idx = start; neighbor_atom_idx < end;
           ++neighbor_atom_idx) {
        if (atom_idx != neighbor_atom_idx && neighbor_atom_idx < end) {
          distance =
            CaculateDistance(
            d_box, shared_px[threadIdx.x], shared_py[threadIdx.x],
            shared_pz[threadIdx.x], __ldg(&px[neighbor_atom_idx]),
            __ldg(&py[neighbor_atom_idx]), __ldg(&pz[neighbor_atom_idx]));
          if (distance < trunc_distance_power_2 &&
              atom_idx < neighbor_atom_idx) {
            neighbors[neighbor_num] = neighbor_atom_idx;
            ++neighbor_num;
          } else if (distance >= trunc_distance_power_2 &&
                     distance < cutoff_2) {
            if (index_shell % selection_frequency == 0 &&
                random_neighbor_num <
                    neighbor_sample_num) {  // random_neighbor_num <
                                            // neighbor_sample_num
              random_neighbors[random_neighbor_start] = neighbor_atom_idx;
              random_neighbor_start++;
              random_neighbor_num++;
            }
            index_shell++;
          }
        }
      }
    }
    neighbor_end[atom_idx] = neighbor_num;
    random_neighbors_num[atom_idx] = random_neighbor_num;
    rbmd::Id my_total_neighbor_num = neighbor_num - neighbor_start[atom_idx];
    if (my_total_neighbor_num > max_neighbor_num[atom_idx]) {
      atomicOr(should_realloc,RBMD_TRUE);
    }
  }
}

void GenerateRblHalfNeighborListOp<device::DEVICE_GPU>::operator()(
    rbmd::Id* per_atom_cell_id, rbmd::Id* in_atom_list_start_index,
    rbmd::Id* in_atom_list_end_index, rbmd::Real trunc_distance_power_2,
    rbmd::Real cutoff_2, rbmd::Id total_atom_num, rbmd::Real* px,
    rbmd::Real* py, rbmd::Real* pz, rbmd::Id* max_neighbor_num,
    rbmd::Id* neighbor_start, rbmd::Id* neighbor_end, rbmd::Id* neighbors,
    rbmd::Id neighbor_sample_num, rbmd::Id* random_neighbors,
    rbmd::Id* random_neighbors_num, Box* d_box, rbmd::Id* should_realloc,
    rbmd::Id* neighbor_cell, rbmd::Id neighbor_cell_num,
    rbmd::Id selection_frequency) {
  unsigned int blocks_per_grid = (total_atom_num + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(
      GenerateRBLHalfNeighborList<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
          per_atom_cell_id, in_atom_list_start_index, in_atom_list_end_index,
          trunc_distance_power_2, cutoff_2, total_atom_num, px, py, pz,
          max_neighbor_num, neighbor_start, neighbor_end, neighbors,
          neighbor_sample_num, random_neighbors, random_neighbors_num, d_box,
          should_realloc, neighbor_cell, neighbor_cell_num,
          selection_frequency));
}
}  // namespace op
