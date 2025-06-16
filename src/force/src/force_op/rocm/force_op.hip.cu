#include "../common/rbmd_define.h"
#include "force/src/force_op/force_op.h"
#include "force_op.h"
#include "model/box.h"

namespace op {

// RBL: Fix RBL Force Single
__global__ void FixRBLForceSingle(const rbmd::Id num_atoms,
                            const rbmd::Real corr_value_x,
                            rbmd::Real* __restrict__ fx) {
  //__shared__ rbmd::Real s_fx[BLOCK_SIZE];

  unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
 // unsigned int local_tid = threadIdx.x;

  //
  if (tid < num_atoms) {
   // s_fx[local_tid] = fx[tid];
   fx[tid] -= corr_value_x;
  }

  // __syncthreads();
  //
  // //
  // if (tid < num_atoms) {
  //   s_fx[local_tid] -= corr_value_x;
  // }
  //
  // __syncthreads();
  //
  // //
  // if (tid < num_atoms) {
  //   fx[tid] = s_fx[local_tid];
  // }
}



// RBL: Fix RBL Force
// __global__ void FixRBLForce(const rbmd::Id num_atoms,
//                             const rbmd::Real corr_value_x,
//                             const rbmd::Real corr_value_y,
//                             const rbmd::Real corr_value_z, rbmd::Real* __restrict__ fx,
//                             rbmd::Real* __restrict__ fy, rbmd::Real* __restrict__ fz) {
//   __shared__ rbmd::Real s_fx[BLOCK_SIZE];
//   __shared__ rbmd::Real s_fy[BLOCK_SIZE];
//   __shared__ rbmd::Real s_fz[BLOCK_SIZE];
//
//   unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
//   unsigned int local_tid = threadIdx.x;
//
//   //
//   if (tid < num_atoms) {
//     s_fx[local_tid] = fx[tid];
//     s_fy[local_tid] = fy[tid];
//     s_fz[local_tid] = fz[tid];
//   }
//
//   __syncthreads();
//
//   //
//   if (tid < num_atoms) {
//     s_fx[local_tid] -= corr_value_x;
//     s_fy[local_tid] -= corr_value_y;
//     s_fz[local_tid] -= corr_value_z;
//   }
//
//   __syncthreads();
//
//   //
//   if (tid < num_atoms) {
//     fx[tid] = s_fx[local_tid];
//     fy[tid] = s_fy[local_tid];
//     fz[tid] = s_fz[local_tid];
//   }
// }



// RBL: Fix LJForce
void FixRBLForceSingleOp<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms, const rbmd::Real corr_value_x,
    rbmd::Real* fx) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

  CHECK_KERNEL(FixRBLForceSingle<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms, corr_value_x,  fx));
}



// RBL: Fix LJForce
// void FixRBLForceOp<device::DEVICE_GPU>::operator()(
//     const rbmd::Id num_atoms, const rbmd::Real corr_value_x,
//     const rbmd::Real corr_value_y, const rbmd::Real corr_value_z,
//     rbmd::Real* fx, rbmd::Real* fy, rbmd::Real* fz) {
//   unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
//
//   CHECK_KERNEL(FixRBLForce<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
//       num_atoms, corr_value_x, corr_value_y, corr_value_z, fx, fy, fz));
// }


}