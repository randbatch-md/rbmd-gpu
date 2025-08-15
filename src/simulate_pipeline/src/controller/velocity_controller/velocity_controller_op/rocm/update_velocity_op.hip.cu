#include "rbmd_define.h"
#include "update_velocity_op.h"
namespace op {
#define THREADS_PER_BLOCK 256

__global__ void UpdateVelocity(const rbmd::Id num_atoms, const rbmd::Real dt,
                               const rbmd::Real fmt2v,
                               const rbmd::Id* atoms_type,
                               const rbmd::Real* mass, const rbmd::Real* fx,
                               const rbmd::Real* fy, const rbmd::Real* fz,
                               rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    rbmd::Id typei = atoms_type[tid];
    rbmd::Real sum_vx = vx[tid];
    rbmd::Real sum_vy = vy[tid];
    rbmd::Real sum_vz = vz[tid];

    sum_vx += 0.5 * fx[tid] / mass[typei] * dt * fmt2v;
    sum_vy += 0.5 * fy[tid] / mass[typei] * dt * fmt2v;
    sum_vz += 0.5 * fz[tid] / mass[typei] * dt * fmt2v;

    vx[tid] = sum_vx;
    vy[tid] = sum_vy;
    vz[tid] = sum_vz;
  }
}

void UpdateVelocityOp<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms, const rbmd::Real dt, const rbmd::Real fmt2v,
    const rbmd::Id* atoms_type, const rbmd::Real* mass, const rbmd::Real* fx,
    const rbmd::Real* fy, const rbmd::Real* fz, rbmd::Real* vx, rbmd::Real* vy,
    rbmd::Real* vz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdateVelocity<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms, dt, fmt2v, atoms_type, mass, fx, fy, fz, vx, vy, vz));
}

__global__ void UpdateVelocity1(const rbmd::Id num_atoms, const rbmd::Real c1,const rbmd::Real dt,
                               const rbmd::Real fmt2v,
                               const rbmd::Id* atoms_type,
                               const rbmd::Real* mass, const rbmd::Real* fx,
                               const rbmd::Real* fy, const rbmd::Real* fz,
                               rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    rbmd::Id typei = atoms_type[tid];
    rbmd::Real sum_vx = vx[tid];
    rbmd::Real sum_vy = vy[tid];
    rbmd::Real sum_vz = vz[tid];

    sum_vx += c1 * fx[tid] / mass[typei] * dt * fmt2v;
    sum_vy += c1 * fy[tid] / mass[typei] * dt * fmt2v;
    sum_vz += c1 * fz[tid] / mass[typei] * dt * fmt2v;

    vx[tid] = sum_vx;
    vy[tid] = sum_vy;
    vz[tid] = sum_vz;
  }
}

void UpdateVelocityOp1<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms, const rbmd::Real c1, const rbmd::Real dt, const rbmd::Real fmt2v,
    const rbmd::Id* atoms_type, const rbmd::Real* mass, const rbmd::Real* fx,
    const rbmd::Real* fy, const rbmd::Real* fz, rbmd::Real* vx, rbmd::Real* vy,
    rbmd::Real* vz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdateVelocity1<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms, c1, dt, fmt2v, atoms_type, mass, fx, fy, fz, vx, vy, vz));
}

__global__ void UpdateVelocity2(const rbmd::Id num_atoms, const rbmd::Real c2, const rbmd::Real dt,
                               const rbmd::Real fmt2v,
                               const rbmd::Id* atoms_type,
                               const rbmd::Real* mass, const rbmd::Real* fx,
                               const rbmd::Real* fy, const rbmd::Real* fz,
                               rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    rbmd::Id typei = atoms_type[tid];
    rbmd::Real sum_vx = vx[tid];
    rbmd::Real sum_vy = vy[tid];
    rbmd::Real sum_vz = vz[tid];

    sum_vx += c2 * fx[tid] / mass[typei] * dt * fmt2v;
    sum_vy += c2 * fy[tid] / mass[typei] * dt * fmt2v;
    sum_vz += c2 * fz[tid] / mass[typei] * dt * fmt2v;

    vx[tid] = sum_vx;
    vy[tid] = sum_vy;
    vz[tid] = sum_vz;
  }
}

void UpdateVelocityOp2<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms, const rbmd::Real c2, const rbmd::Real dt, const rbmd::Real fmt2v,
    const rbmd::Id* atoms_type, const rbmd::Real* mass, const rbmd::Real* fx,
    const rbmd::Real* fy, const rbmd::Real* fz, rbmd::Real* vx, rbmd::Real* vy,
    rbmd::Real* vz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdateVelocity2<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms,c2, dt, fmt2v, atoms_type, mass, fx, fy, fz, vx, vy, vz));
}

__global__ void UpdateVelocity3(const rbmd::Id num_atoms, const rbmd::Real c3, const rbmd::Real dt,
                               const rbmd::Real fmt2v,
                               const rbmd::Id* atoms_type,
                               const rbmd::Real* mass, const rbmd::Real* fx,
                               const rbmd::Real* fy, const rbmd::Real* fz,
                               rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    rbmd::Id typei = atoms_type[tid];
    rbmd::Real sum_vx = vx[tid];
    rbmd::Real sum_vy = vy[tid];
    rbmd::Real sum_vz = vz[tid];

    sum_vx += c3 * fx[tid] / mass[typei] * dt * fmt2v;
    sum_vy += c3 * fy[tid] / mass[typei] * dt * fmt2v;
    sum_vz += c3 * fz[tid] / mass[typei] * dt * fmt2v;

    vx[tid] = sum_vx;
    vy[tid] = sum_vy;
    vz[tid] = sum_vz;
  }
}

void UpdateVelocityOp3<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms, const rbmd::Real c3, const rbmd::Real dt, const rbmd::Real fmt2v,
    const rbmd::Id* atoms_type, const rbmd::Real* mass, const rbmd::Real* fx,
    const rbmd::Real* fy, const rbmd::Real* fz, rbmd::Real* vx, rbmd::Real* vy,
    rbmd::Real* vz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdateVelocity3<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms,c3, dt, fmt2v, atoms_type, mass, fx, fy, fz, vx, vy, vz));
}

__global__ void UpdateVelocity4(const rbmd::Id num_atoms, const rbmd::Real c4, const rbmd::Real dt,
                               const rbmd::Real fmt2v,
                               const rbmd::Id* atoms_type,
                               const rbmd::Real* mass, const rbmd::Real* fx,
                               const rbmd::Real* fy, const rbmd::Real* fz,
                               rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    rbmd::Id typei = atoms_type[tid];
    rbmd::Real sum_vx = vx[tid];
    rbmd::Real sum_vy = vy[tid];
    rbmd::Real sum_vz = vz[tid];

    sum_vx += c4 * fx[tid] / mass[typei] * dt * fmt2v;
    sum_vy += c4 * fy[tid] / mass[typei] * dt * fmt2v;
    sum_vz += c4 * fz[tid] / mass[typei] * dt * fmt2v;

    vx[tid] = sum_vx;
    vy[tid] = sum_vy;
    vz[tid] = sum_vz;
  }
}

void UpdateVelocityOp4<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms, const rbmd::Real c4, const rbmd::Real dt, const rbmd::Real fmt2v,
    const rbmd::Id* atoms_type, const rbmd::Real* mass, const rbmd::Real* fx,
    const rbmd::Real* fy, const rbmd::Real* fz, rbmd::Real* vx, rbmd::Real* vy,
    rbmd::Real* vz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdateVelocity4<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms,c4, dt, fmt2v, atoms_type, mass, fx, fy, fz, vx, vy, vz));
}
}  // namespace op
