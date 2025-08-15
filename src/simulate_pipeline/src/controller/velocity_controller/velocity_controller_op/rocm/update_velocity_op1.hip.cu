#include "rbmd_define.h"
#include "update_velocity_op.h"
namespace op {
#define THREADS_PER_BLOCK 256

__global__ void UpdateVelocityvl(const rbmd::Id num_atoms, const rbmd::Real dt,
                               const rbmd::Real fmt2v,
                               const rbmd::Id* atoms_type,
                               const rbmd::Real* mass, const rbmd::Real* fx,
                               const rbmd::Real* fy, const rbmd::Real* fz,
                               const rbmd::Real* px, const rbmd::Real* py, const rbmd::Real* pz,
                               rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    rbmd::Id typei = atoms_type[tid];
    rbmd::Real sum_vx = vx[tid];
    rbmd::Real sum_vy = vy[tid];
    rbmd::Real sum_vz = vz[tid];
    rbmd::Real cont = 0.50000;

    sum_vx +=  cont * fx[tid] / mass[typei] * dt * fmt2v;
    sum_vy +=  cont * fy[tid] / mass[typei] * dt * fmt2v;
    sum_vz +=  cont * fz[tid] / mass[typei] * dt * fmt2v;

    vx[tid] = sum_vx;
    vy[tid] = sum_vy;
    vz[tid] = sum_vz;

  }
}

void UpdateVelocityOpvl<device::DEVICE_GPU>::operator()(
const rbmd::Id num_atoms, const rbmd::Real dt,
                          const rbmd::Real fmt2v,
                          const rbmd::Id* atoms_type,
                          const rbmd::Real* mass, const rbmd::Real* fx,
                          const rbmd::Real* fy, const rbmd::Real* fz,
                          const rbmd::Real* px, const rbmd::Real* py, const rbmd::Real* pz,
                          rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdateVelocityvl<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms, dt, fmt2v, atoms_type, mass, fx, fy, fz, px, py, pz, vx, vy, vz));
}

__global__ void UpdateVelocitybm(const rbmd::Id num_atoms, const rbmd::Real dt,
  rbmd::Id test_current_step,const rbmd::Real fmt2v,const rbmd::Id* atoms_type,
  const rbmd::Real* mass, const rbmd::Real* fx,const rbmd::Real* fy, const rbmd::Real* fz,
  rbmd::Real* prev_fx, rbmd::Real* prev_fy, rbmd::Real* prev_fz,
  rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {

    if (test_current_step  < 2){
      rbmd::Id typei = atoms_type[tid];
      rbmd::Real sum_vx = vx[tid];
      rbmd::Real sum_vy = vy[tid];
      rbmd::Real sum_vz = vz[tid];
      rbmd::Real cont = 0.50000;

      sum_vx +=  cont * fx[tid] / mass[typei] * dt * fmt2v;
      sum_vy +=  cont * fy[tid] / mass[typei] * dt * fmt2v;
      sum_vz +=  cont * fz[tid] / mass[typei] * dt * fmt2v;

      vx[tid] = sum_vx;
      vy[tid] = sum_vy;
      vz[tid] = sum_vz;

      prev_fx[tid] = fx[tid];
      prev_fy[tid] = fy[tid];
      prev_fz[tid] = fz[tid];
    }
    else {
      rbmd::Id typei = atoms_type[tid];
      rbmd::Real sum_vx = vx[tid];
      rbmd::Real sum_vy = vy[tid];
      rbmd::Real sum_vz = vz[tid];
      rbmd::Real cont3 = 0.16666667;

      sum_vx += cont3 * (2 * fx[tid] + prev_fx[tid]) / mass[typei] * dt * fmt2v;
      sum_vy += cont3 * (2 * fy[tid] + prev_fy[tid]) / mass[typei] * dt * fmt2v;
      sum_vz += cont3 * (2 * fz[tid] + prev_fz[tid]) / mass[typei] * dt * fmt2v;

      vx[tid] = sum_vx;
      vy[tid] = sum_vy;
      vz[tid] = sum_vz;

      prev_fx[tid] = fx[tid];
      prev_fy[tid] = fy[tid];
      prev_fz[tid] = fz[tid];
    }
  }
}

void UpdateVelocityOpbm<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms, const rbmd::Real dt, rbmd::Id test_current_step,
    const rbmd::Real fmt2v,const rbmd::Id* atoms_type,const rbmd::Real* mass,
    const rbmd::Real* fx,const rbmd::Real* fy, const rbmd::Real* fz,
    rbmd::Real* prev_fx, rbmd::Real* prev_fy, rbmd::Real* prev_fz,
    rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdateVelocitybm<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms, dt, test_current_step, fmt2v, atoms_type, mass, fx, fy, fz, prev_fx, prev_fy, prev_fz, vx, vy, vz));
}

}  // namespace op
