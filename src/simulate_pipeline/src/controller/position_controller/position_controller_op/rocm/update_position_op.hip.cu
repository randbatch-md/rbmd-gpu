#include "rbmd_define.h"
#include "update_position_op.h"

namespace op {
#define THREADS_PER_BLOCK 256

__global__ void UpdatePositionFlag0(
    const rbmd::Id num_atoms, const rbmd::Real dt, Box  box  ,
    const rbmd::Real* vx, const rbmd::Real* vy, const rbmd::Real* vz,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, rbmd::Id* flag_px,
    rbmd::Id* flag_py, rbmd::Id* flag_pz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    rbmd::Real sum_px = px[tid];
    rbmd::Real sum_py = py[tid];
    rbmd::Real sum_pz = pz[tid];

    sum_px += vx[tid] * dt;
    sum_py += vy[tid] * dt;
    sum_pz += vz[tid] * dt;

    px[tid] = sum_px;
    py[tid] = sum_py;
    pz[tid] = sum_pz;

    ApplyPBC(box, px[tid], py[tid], pz[tid],
      flag_px[tid], flag_py[tid],flag_pz[tid]);
  }
}

__global__ void UpdatePositionFlag(
    const rbmd::Id num_atoms, const rbmd::Real dt, Box box,
    const rbmd::Real* vx, const rbmd::Real* vy, const rbmd::Real* vz,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz,
    rbmd::Id* flag_px, rbmd::Id* flag_py, rbmd::Id* flag_pz)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    // 1. 计算不受约束的、更新后的“真实”坐标
    rbmd::Real new_px = px[tid] + vx[tid] * dt;
    rbmd::Real new_py = py[tid] + vy[tid] * dt;
    rbmd::Real new_pz = pz[tid] + vz[tid] * dt;

    // 2. 读取当前的映像标志
    rbmd::Id current_flag_x = flag_px[tid];
    rbmd::Id current_flag_y = flag_py[tid];
    rbmd::Id current_flag_z = flag_pz[tid];

    // 3. 调用鲁棒的PBC函数，它会同时更新坐标和映像标志
    ApplyPBC_Robust(box, new_px, new_py, new_pz,
                    current_flag_x, current_flag_y, current_flag_z);

    // 4. 将更新后的值写回全局内存
    px[tid] = new_px;
    py[tid] = new_py;
    pz[tid] = new_pz;
    flag_px[tid] = current_flag_x;
    flag_py[tid] = current_flag_y;
    flag_pz[tid] = current_flag_z;
  }
}

__global__ void UpdatePosition(const rbmd::Id num_atoms, const rbmd::Real dt, Box  box  ,
                               const rbmd::Real* vx, const rbmd::Real* vy,
                               const rbmd::Real* vz, rbmd::Real* px,
                               rbmd::Real* py, rbmd::Real* pz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    rbmd::Real sum_px = px[tid];
    rbmd::Real sum_py = py[tid];
    rbmd::Real sum_pz = pz[tid];

    sum_px += vx[tid] * dt;
    sum_py += vy[tid] * dt;
    sum_pz += vz[tid] * dt;

    px[tid] = sum_px;
    py[tid] = sum_py;
    pz[tid] = sum_pz;

    ApplyPBC_unflag(box, px[tid], py[tid], pz[tid]);
  }
}

__global__ void PBC(
    const rbmd::Id num_atoms, Box  box  ,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, rbmd::Id* flag_px,
    rbmd::Id* flag_py, rbmd::Id* flag_pz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    ApplyPBC_Robust(box, px[tid], py[tid], pz[tid],
      flag_px[tid], flag_py[tid],flag_pz[tid]);
  }
}


void UpdatePositionFlagOp<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms, const rbmd::Real dt, Box  box  ,
    const rbmd::Real* vx, const rbmd::Real* vy, const rbmd::Real* vz,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, rbmd::Id* flag_px,
    rbmd::Id* flag_py, rbmd::Id* flag_pz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdatePositionFlag<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms, dt, box, vx, vy, vz, px, py, pz, flag_px, flag_py, flag_pz));
}

void UpdatePositionOp<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms, const rbmd::Real dt, Box  box  , const rbmd::Real* vx,
    const rbmd::Real* vy, const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
    rbmd::Real* pz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdatePosition<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms, dt, box, vx, vy, vz, px,py, pz));
}

void PBCOp<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms,  Box  box  ,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, rbmd::Id* flag_px,
    rbmd::Id* flag_py, rbmd::Id* flag_pz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(PBC<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms, box, px, py, pz, flag_px, flag_py, flag_pz));
}

}  // namespace op
