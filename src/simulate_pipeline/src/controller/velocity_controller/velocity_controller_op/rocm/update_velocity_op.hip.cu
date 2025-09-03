#include "rbmd_define.h"
#include "update_velocity_op.h"
namespace op {
#define THREADS_PER_BLOCK 256

// 蛙跳法leapfrog
__global__ void UpdateVelocity(const rbmd::Id num_atoms, const rbmd::Real dt,
                               const rbmd::Real fmt2v,
                               const rbmd::Id* atoms_type,
                               const rbmd::Real* mass, const rbmd::Real* fx,
                               const rbmd::Real* fy, const rbmd::Real* fz,
                               rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    // 获取当前原子的类型ID
    rbmd::Id typei = atoms_type[tid];
    rbmd::Real sum_vx = vx[tid];
    rbmd::Real sum_vy = vy[tid];
    rbmd::Real sum_vz = vz[tid];

    // 更新速度
    sum_vx += 0.5 * fx[tid] / mass[typei] * fmt2v * dt;
    sum_vy += 0.5 * fy[tid] / mass[typei] * fmt2v * dt;
    sum_vz += 0.5 * fz[tid] / mass[typei] * fmt2v * dt;

    // 更新速度写回全局内存
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

// vv
__global__ void UpdateVelocityvv(const rbmd::Id num_atoms, const rbmd::Real dt,
                               const rbmd::Real fmt2v, const rbmd::Id* atoms_type, const rbmd::Real* mass,
                               const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz,
                               const rbmd::Real* fx_prev, const rbmd::Real* fy_prev, const rbmd::Real* fz_prev,
                               rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    // 获取当前原子的类型ID
    rbmd::Id typei = atoms_type[tid];
    rbmd::Real sum_vx = vx[tid];
    rbmd::Real sum_vy = vy[tid];
    rbmd::Real sum_vz = vz[tid];
    rbmd::Real const_half = 0.5000;

    // 更新速度
    sum_vx += const_half * (fx[tid] + fx_prev[tid]) / mass[typei] * fmt2v * dt;
    sum_vy += const_half * (fy[tid] + fy_prev[tid]) / mass[typei] * fmt2v * dt;
    sum_vz += const_half * (fz[tid] + fz_prev[tid]) / mass[typei] * fmt2v * dt;

    // 更新速度写回全局内存
    vx[tid] = sum_vx;
    vy[tid] = sum_vy;
    vz[tid] = sum_vz;
  }
}

void UpdateVelocityOpvv<device::DEVICE_GPU>::operator()(
                           const rbmd::Id num_atoms, const rbmd::Real dt,
                           const rbmd::Real fmt2v, const rbmd::Id* atoms_type, const rbmd::Real* mass,
                           const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz,
                           const rbmd::Real* fx_prev, const rbmd::Real* fy_prev, const rbmd::Real* fz_prev,
                           rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdateVelocityvv<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms, dt, fmt2v, atoms_type, mass, fx, fy, fz,fx_prev,fy_prev,fz_prev, vx, vy, vz));
}

// PRK
__global__ void UpdateVelocity1(const rbmd::Id num_atoms, const rbmd::Real c1, const rbmd::Real dt,
                               const rbmd::Real fmt2v,
                               const rbmd::Id* atoms_type,
                               const rbmd::Real* mass, const rbmd::Real* fx,
                               const rbmd::Real* fy, const rbmd::Real* fz,
                               rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;


  if (tid < num_atoms) {
    // 获取当前原子的类型ID
    rbmd::Id typei = atoms_type[tid];
    rbmd::Real sum_vx = vx[tid];
    rbmd::Real sum_vy = vy[tid];
    rbmd::Real sum_vz = vz[tid];


    // 更新速度
    sum_vx += c1 * fx[tid] / mass[typei] * fmt2v * dt;
    sum_vy += c1 * fy[tid] / mass[typei] * fmt2v * dt;
    sum_vz += c1 * fz[tid] / mass[typei] * fmt2v * dt;

    // 更新速度写回全局内存
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
    // 获取当前原子的类型ID
    rbmd::Id typei = atoms_type[tid];
    rbmd::Real sum_vx = vx[tid];
    rbmd::Real sum_vy = vy[tid];
    rbmd::Real sum_vz = vz[tid];

    // 更新速度
    sum_vx += c2 * fx[tid] / mass[typei] * fmt2v * dt;
    sum_vy += c2 * fy[tid] / mass[typei] * fmt2v * dt;
    sum_vz += c2 * fz[tid] / mass[typei] * fmt2v * dt;

    // 更新速度写回全局内存
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
      num_atoms, c2, dt, fmt2v, atoms_type, mass, fx, fy, fz, vx, vy, vz));
}


__global__ void UpdateVelocity3(const rbmd::Id num_atoms, const rbmd::Real c3, const rbmd::Real dt,
                               const rbmd::Real fmt2v,
                               const rbmd::Id* atoms_type,
                               const rbmd::Real* mass, const rbmd::Real* fx,
                               const rbmd::Real* fy, const rbmd::Real* fz,
                               rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    // 获取当前原子的类型ID
    rbmd::Id typei = atoms_type[tid];
    rbmd::Real sum_vx = vx[tid];
    rbmd::Real sum_vy = vy[tid];
    rbmd::Real sum_vz = vz[tid];

    // 更新速度
    sum_vx += c3 * fx[tid] / mass[typei] * fmt2v * dt;
    sum_vy += c3 * fy[tid] / mass[typei] * fmt2v * dt;
    sum_vz += c3 * fz[tid] / mass[typei] * fmt2v * dt;

    // 更新速度写回全局内存
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
      num_atoms, c3, dt, fmt2v, atoms_type, mass, fx, fy, fz, vx, vy, vz));
}


__global__ void UpdateVelocity4(const rbmd::Id num_atoms, const rbmd::Real c4, const rbmd::Real dt,
                               const rbmd::Real fmt2v,
                               const rbmd::Id* atoms_type,
                               const rbmd::Real* mass, const rbmd::Real* fx,
                               const rbmd::Real* fy, const rbmd::Real* fz,
                               rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    // 获取当前原子的类型ID
    rbmd::Id typei = atoms_type[tid];
    rbmd::Real sum_vx = vx[tid];
    rbmd::Real sum_vy = vy[tid];
    rbmd::Real sum_vz = vz[tid];

    // 更新速度
    sum_vx += c4 * fx[tid] / mass[typei] * fmt2v * dt;
    sum_vy += c4 * fy[tid] / mass[typei] * fmt2v * dt;
    sum_vz += c4 * fz[tid] / mass[typei] * fmt2v * dt;

    // 更新速度写回全局内存
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
      num_atoms, c4, dt, fmt2v, atoms_type, mass, fx, fy, fz, vx, vy, vz));
}

// Beeman
__global__ void UpdateVelocityBeeman(
    const rbmd::Id num_atoms, const rbmd::Real dt, rbmd::Id test_current_step,const rbmd::Real fmt2v,
    const rbmd::Id* atoms_type, const rbmd::Real* mass,
    const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz,
    rbmd::Real* f_pre1_x, rbmd::Real* f_pre1_y, rbmd::Real* f_pre1_z,
    rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid < num_atoms) {
    if (test_current_step < 2){
      rbmd::Id typei = atoms_type[tid];
      rbmd::Real sum_vx = vx[tid];
      rbmd::Real sum_vy = vy[tid];
      rbmd::Real sum_vz = vz[tid];
      rbmd::Real const_half = 0.5000;

      sum_vx +=  const_half * fx[tid] / mass[typei] * dt * fmt2v;
      sum_vy +=  const_half * fy[tid] / mass[typei] * dt * fmt2v;
      sum_vz +=  const_half * fz[tid] / mass[typei] * dt * fmt2v;

      f_pre1_x[tid] = fx[tid];
      f_pre1_y[tid] = fy[tid];
      f_pre1_z[tid] = fz[tid];

      vx[tid] = sum_vx;
      vy[tid] = sum_vy;
      vz[tid] = sum_vz;
    }
    else {
      rbmd::Id typei = atoms_type[tid];
      rbmd::Real sum_vx = vx[tid];
      rbmd::Real sum_vy = vy[tid];
      rbmd::Real sum_vz = vz[tid];
      rbmd::Real const_one_third = 0.33333333;

      sum_vx += const_one_third * (fx[tid]) / mass[typei] * dt * fmt2v;
      sum_vy += const_one_third * (fy[tid]) / mass[typei] * dt * fmt2v;
      sum_vz += const_one_third * (fz[tid]) / mass[typei] * dt * fmt2v;

      // sum_vx += const_one_third * (fx[tid] + 0.5*f_pre1_x[tid]) / mass[typei] * dt * fmt2v;
      // sum_vy += const_one_third * (fy[tid] + 0.5*f_pre1_y[tid]) / mass[typei] * dt * fmt2v;
      // sum_vz += const_one_third * (fz[tid] + 0.5*f_pre1_z[tid]) / mass[typei] * dt * fmt2v;

      vx[tid] = sum_vx;
      vy[tid] = sum_vy;
      vz[tid] = sum_vz;

      f_pre1_x[tid] = fx[tid];
      f_pre1_y[tid] = fy[tid];
      f_pre1_z[tid] = fz[tid];
    }
  }
}
void UpdateVelocityOpBeeman<device::DEVICE_GPU>::operator()(
  const rbmd::Id num_atoms, const rbmd::Real dt, rbmd::Id test_current_step, const rbmd::Real fmt2v,
  const rbmd::Id* atoms_type, const rbmd::Real* mass,
  const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz,                 // F(t+Δt)
  rbmd::Real* f_pre1_x, rbmd::Real* f_pre1_y, rbmd::Real* f_pre1_z, // F(t-Δt)
  rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdateVelocityBeeman<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms, dt, test_current_step, fmt2v, atoms_type, mass,
      fx, fy, fz,
      f_pre1_x, f_pre1_y, f_pre1_z,
      vx, vy, vz));
}

}// namespace op

