#include "rbmd_define.h"
#include "update_position_op.h"

namespace op {
#define THREADS_PER_BLOCK 256

__global__ void UpdatePositionFlag(
    const rbmd::Id num_atoms, const rbmd::Real dt,const rbmd::Real fmt2v, Box  box  , const rbmd::Id* atoms_type,
    const rbmd::Real* vx, const rbmd::Real* vy, const rbmd::Real* vz, const rbmd::Real* mass,
    const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, rbmd::Id* flag_px,
    rbmd::Id* flag_py, rbmd::Id* flag_pz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    // rbmd::Id typei = atoms_type[tid];
    rbmd::Real sum_px = px[tid];
    rbmd::Real sum_py = py[tid];
    rbmd::Real sum_pz = pz[tid];

    sum_px += vx[tid] * dt; // + 0.5 * fx[tid] / mass[typei] * dt * dt * fmt2v;
    sum_py += vy[tid] * dt; // + 0.5 * fy[tid] / mass[typei] * dt * dt * fmt2v;
    sum_pz += vz[tid] * dt; // + 0.5 * fz[tid] / mass[typei] * dt * dt * fmt2v;

    px[tid] = sum_px;
    py[tid] = sum_py;
    pz[tid] = sum_pz;

    ApplyPBC(box, px[tid], py[tid], pz[tid],
      flag_px[tid], flag_py[tid],flag_pz[tid]);
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

void UpdatePositionFlagOp<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms, const rbmd::Real dt, const rbmd::Real fmt2v, Box  box , const rbmd::Id* atoms_type,
    const rbmd::Real* vx, const rbmd::Real* vy, const rbmd::Real* vz, const rbmd::Real* mass,
    const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, rbmd::Id* flag_px,
    rbmd::Id* flag_py, rbmd::Id* flag_pz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdatePositionFlag<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms, dt, fmt2v, box, atoms_type, vx, vy, vz, mass, fx, fy, fz, px, py, pz, flag_px, flag_py, flag_pz));
}

// void UpdatePositionOp<device::DEVICE_GPU>::operator()(
//     const rbmd::Id num_atoms, const rbmd::Real dt, Box  box  , const rbmd::Real* vx,
//     const rbmd::Real* vy, const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
//     rbmd::Real* pz) {
//   unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
//   CHECK_KERNEL(UpdatePosition<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
//       num_atoms, dt, box, vx, vy, vz, px,py, pz));
// }

__global__ void UpdatePositionFlag1(
    const rbmd::Id num_atoms,const rbmd::Real d1, const rbmd::Real dt, Box  box  ,
    const rbmd::Real* vx, const rbmd::Real* vy, const rbmd::Real* vz,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, rbmd::Id* flag_px,
    rbmd::Id* flag_py, rbmd::Id* flag_pz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    rbmd::Real sum_px = px[tid];
    rbmd::Real sum_py = py[tid];
    rbmd::Real sum_pz = pz[tid];

    sum_px += d1 * vx[tid] * dt;
    sum_py += d1 * vy[tid] * dt;
    sum_pz += d1 * vz[tid] * dt;

    px[tid] = sum_px;
    py[tid] = sum_py;
    pz[tid] = sum_pz;

    ApplyPBC(box, px[tid], py[tid], pz[tid],
      flag_px[tid], flag_py[tid],flag_pz[tid]);
  }
}


void UpdatePositionFlagOp1<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms,const rbmd::Real d1, const rbmd::Real dt, Box  box  ,
    const rbmd::Real* vx, const rbmd::Real* vy, const rbmd::Real* vz,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, rbmd::Id* flag_px,
    rbmd::Id* flag_py, rbmd::Id* flag_pz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdatePositionFlag1<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms,d1, dt, box, vx, vy, vz, px, py, pz, flag_px, flag_py, flag_pz));
}
__global__ void UpdatePositionFlag2(
    const rbmd::Id num_atoms,const rbmd::Real d2, const rbmd::Real dt, Box  box  ,
    const rbmd::Real* vx, const rbmd::Real* vy, const rbmd::Real* vz,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, rbmd::Id* flag_px,
    rbmd::Id* flag_py, rbmd::Id* flag_pz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    rbmd::Real sum_px = px[tid];
    rbmd::Real sum_py = py[tid];
    rbmd::Real sum_pz = pz[tid];

    sum_px += d2 * vx[tid] * dt;
    sum_py += d2 * vy[tid] * dt;
    sum_pz += d2 * vz[tid] * dt;

    px[tid] = sum_px;
    py[tid] = sum_py;
    pz[tid] = sum_pz;

    ApplyPBC(box, px[tid], py[tid], pz[tid],
      flag_px[tid], flag_py[tid],flag_pz[tid]);
  }
}


void UpdatePositionFlagOp2<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms,const rbmd::Real d2, const rbmd::Real dt, Box  box  ,
    const rbmd::Real* vx, const rbmd::Real* vy, const rbmd::Real* vz,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, rbmd::Id* flag_px,
    rbmd::Id* flag_py, rbmd::Id* flag_pz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdatePositionFlag2<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms,d2, dt, box, vx, vy, vz, px, py, pz, flag_px, flag_py, flag_pz));
}
__global__ void UpdatePositionFlag3(
    const rbmd::Id num_atoms,const rbmd::Real d3, const rbmd::Real dt, Box  box  ,
    const rbmd::Real* vx, const rbmd::Real* vy, const rbmd::Real* vz,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, rbmd::Id* flag_px,
    rbmd::Id* flag_py, rbmd::Id* flag_pz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    rbmd::Real sum_px = px[tid];
    rbmd::Real sum_py = py[tid];
    rbmd::Real sum_pz = pz[tid];

    sum_px += d3 * vx[tid] * dt;
    sum_py += d3 * vy[tid] * dt;
    sum_pz += d3 * vz[tid] * dt;

    px[tid] = sum_px;
    py[tid] = sum_py;
    pz[tid] = sum_pz;

    ApplyPBC(box, px[tid], py[tid], pz[tid],
      flag_px[tid], flag_py[tid],flag_pz[tid]);
  }
}


void UpdatePositionFlagOp3<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms,const rbmd::Real d3, const rbmd::Real dt, Box  box  ,
    const rbmd::Real* vx, const rbmd::Real* vy, const rbmd::Real* vz,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, rbmd::Id* flag_px,
    rbmd::Id* flag_py, rbmd::Id* flag_pz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdatePositionFlag3<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms,d3, dt, box, vx, vy, vz, px, py, pz, flag_px, flag_py, flag_pz));
}
__global__ void UpdatePositionFlag4(
    const rbmd::Id num_atoms,const rbmd::Real d4, const rbmd::Real dt, Box  box  ,
    const rbmd::Real* vx, const rbmd::Real* vy, const rbmd::Real* vz,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, rbmd::Id* flag_px,
    rbmd::Id* flag_py, rbmd::Id* flag_pz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    rbmd::Real sum_px = px[tid];
    rbmd::Real sum_py = py[tid];
    rbmd::Real sum_pz = pz[tid];

    sum_px += d4 * vx[tid] * dt;
    sum_py += d4 * vy[tid] * dt;
    sum_pz += d4 * vz[tid] * dt;

    px[tid] = sum_px;
    py[tid] = sum_py;
    pz[tid] = sum_pz;

    ApplyPBC(box, px[tid], py[tid], pz[tid],
      flag_px[tid], flag_py[tid],flag_pz[tid]);
  }
}


void UpdatePositionFlagOp4<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms,const rbmd::Real d4, const rbmd::Real dt, Box  box  ,
    const rbmd::Real* vx, const rbmd::Real* vy, const rbmd::Real* vz,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, rbmd::Id* flag_px,
    rbmd::Id* flag_py, rbmd::Id* flag_pz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdatePositionFlag4<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms,d4, dt, box, vx, vy, vz, px, py, pz, flag_px, flag_py, flag_pz));
}




}  // namespace op
