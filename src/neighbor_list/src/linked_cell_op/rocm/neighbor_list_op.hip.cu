#include <hip/hip_runtime.h>

#include <hipcub/hipcub.hpp>

#include "../../../../data_manager/include/model/box.h"
#include "common/device_types.h"
#include "common/rbmd_define.h"
#include "common/types.h"
#include "neighbor_list_op.h"

namespace op {
__global__ void InitEndIndex(rbmd::Id* neighbor_num, rbmd::Id* start_index,
                             rbmd::Id* end_index, rbmd::Id total_atom_num) {
  const unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < total_atom_num) {
    end_index[idx] = start_index[idx] + neighbor_num[idx];
  }
}

void InitEndIndexOp<device::DEVICE_GPU>::operator()(rbmd::Id* neighbor_num, rbmd::Id* start_index,
                  rbmd::Id* end_index, rbmd::Id total_atom_num) {
    unsigned int blocks_per_grid =
        (total_atom_num + BLOCK_SIZE - 1) / BLOCK_SIZE;
    CHECK_KERNEL(InitEndIndex<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
        neighbor_num, start_index, end_index, total_atom_num));
  }
};

