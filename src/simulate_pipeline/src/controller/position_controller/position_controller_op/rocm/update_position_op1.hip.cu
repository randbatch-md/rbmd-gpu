#include "rbmd_define.h"
#include "update_position_op.h"

namespace op {
#define THREADS_PER_BLOCK 256


__global__ void UpdatePositionFlagvl(
    const rbmd::Id num_atoms,const rbmd::Real fmt2v, const rbmd::Real dt, rbmd::Id test_current_step,Box  box  , const rbmd::Id* atoms_type,
    const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz,const rbmd::Real* mass,
    rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz,
    rbmd::Real* prev_px,rbmd::Real* prev_py,rbmd::Real* prev_pz,
    rbmd::Id* flag_px,rbmd::Id* flag_py, rbmd::Id* flag_pz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    rbmd::Id typei = atoms_type[tid];
    rbmd::Real sum_px = px[tid];
    rbmd::Real sum_py = py[tid];
    rbmd::Real sum_pz = pz[tid];
    rbmd::Real sum_px_prev = prev_px[tid];
    rbmd::Real sum_py_prev = prev_py[tid];
    rbmd::Real sum_pz_prev = prev_pz[tid];
    rbmd::Real sum_vx = vx[tid];
    rbmd::Real sum_vy = vy[tid];
    rbmd::Real sum_vz = vz[tid];
    rbmd::Real cont1 = 0.50000;
    rbmd::Real cont2 = 2.00000;
    if (test_current_step < 1000){

      sum_px += vx[tid] * dt + cont1 * fx[tid] / mass[typei] * dt * dt * fmt2v;
      sum_py += vy[tid] * dt + cont1 * fy[tid] / mass[typei] * dt * dt * fmt2v;
      sum_pz += vz[tid] * dt + cont1 * fz[tid] / mass[typei] * dt * dt * fmt2v;

      sum_vx +=  cont1 * fx[tid] / mass[typei] * dt * fmt2v;
      sum_vy +=  cont1 * fy[tid] / mass[typei] * dt * fmt2v;
      sum_vz +=  cont1 * fz[tid] / mass[typei] * dt * fmt2v;

      prev_px[tid] = px[tid];
      prev_py[tid] = py[tid];
      prev_pz[tid] = pz[tid];

      px[tid] = sum_px;
      py[tid] = sum_py;
      pz[tid] = sum_pz;

      vx[tid] = sum_vx;
      vy[tid] = sum_vy;
      vz[tid] = sum_vz;
    }
    else {
      sum_px = cont2 * sum_px - sum_px_prev + fx[tid] / mass[typei] * dt * dt * fmt2v;
      sum_py = cont2 * sum_px - sum_py_prev + fy[tid] / mass[typei] * dt * dt * fmt2v;
      sum_pz = cont2 * sum_px - sum_pz_prev + fz[tid] / mass[typei] * dt * dt * fmt2v;

      sum_vx +=  cont1 * fx[tid] / mass[typei] * dt * fmt2v;
      sum_vy +=  cont1 * fy[tid] / mass[typei] * dt * fmt2v;
      sum_vz +=  cont1 * fz[tid] / mass[typei] * dt * fmt2v;

      prev_px[tid] = px[tid];
      prev_py[tid] = py[tid];
      prev_pz[tid] = pz[tid];

      px[tid] = sum_px;
      py[tid] = sum_py;
      pz[tid] = sum_pz;

      vx[tid] = sum_vx;
      vy[tid] = sum_vy;
      vz[tid] = sum_vz;
    }
    ApplyPBC(box, px[tid], py[tid], pz[tid],
    flag_px[tid], flag_py[tid],flag_pz[tid]);

  }
}

void UpdatePositionFlagOpvl<device::DEVICE_GPU>::operator()(
   const rbmd::Id num_atoms, const rbmd::Real fmt2v,const rbmd::Real dt, const rbmd::Id test_current_step, Box  box  ,const rbmd::Id* atoms_type,
   const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz,const rbmd::Real* mass,
   rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz,
   rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz,
   rbmd::Real* prev_px,rbmd::Real* prev_py,rbmd::Real* prev_pz,
   rbmd::Id* flag_px,rbmd::Id* flag_py, rbmd::Id* flag_pz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdatePositionFlagvl<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms,fmt2v, dt, test_current_step, box, atoms_type, fx, fy, fz, mass, vx, vy, vz, px, py, pz, prev_px, prev_py, prev_pz, flag_px, flag_py, flag_pz));
}

__global__ void UpdatePositionFlagbm(
    const rbmd::Id num_atoms,const rbmd::Real fmt2v, const rbmd::Real dt, rbmd::Id test_current_step,Box  box  , const rbmd::Id* atoms_type,
    const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz,const rbmd::Real* mass,
    rbmd::Real* prev_fx, rbmd::Real* prev_fy, rbmd::Real* prev_fz,
    rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz,
    rbmd::Id* flag_px,rbmd::Id* flag_py, rbmd::Id* flag_pz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    if (test_current_step < 2){
      rbmd::Id typei = atoms_type[tid];
      rbmd::Real sum_px = px[tid];
      rbmd::Real sum_py = py[tid];
      rbmd::Real sum_pz = pz[tid];
      rbmd::Real sum_vx = vx[tid];
      rbmd::Real sum_vy = vy[tid];
      rbmd::Real sum_vz = vz[tid];
      rbmd::Real cont1 = 0.50000;

      sum_px += vx[tid] * dt + cont1 * fx[tid] / mass[typei] * dt * dt * fmt2v;
      sum_py += vy[tid] * dt + cont1 * fy[tid] / mass[typei] * dt * dt * fmt2v;
      sum_pz += vz[tid] * dt + cont1 * fz[tid] / mass[typei] * dt * dt * fmt2v;

      // if (tid == 1) {
      //     printf("fx[1] = %f, fx_prev[1] = %f\n", fx[1], prev_fx[1]);
      // }
      prev_fx[tid] = fx[tid];
      prev_fy[tid] = fy[tid];
      prev_fz[tid] = fz[tid];
      // if (tid == 1) {
      //   printf("fx_prev[1] = %f\n", prev_fx[1]);
      // }

      sum_vx +=  cont1 * prev_fx[tid] / mass[typei] * dt * fmt2v;
      sum_vy +=  cont1 * prev_fy[tid] / mass[typei] * dt * fmt2v;
      sum_vz +=  cont1 * prev_fz[tid] / mass[typei] * dt * fmt2v;


      px[tid] = sum_px;
      py[tid] = sum_py;
      pz[tid] = sum_pz;

      vx[tid] = sum_vx;
      vy[tid] = sum_vy;
      vz[tid] = sum_vz;

      ApplyPBC(box, px[tid], py[tid], pz[tid],
        flag_px[tid], flag_py[tid],flag_pz[tid]);
    }
    else {
      rbmd::Id typei = atoms_type[tid];
      rbmd::Real sum_px = px[tid];
      rbmd::Real sum_py = py[tid];
      rbmd::Real sum_pz = pz[tid];
      rbmd::Real sum_vx = vx[tid];
      rbmd::Real sum_vy = vy[tid];
      rbmd::Real sum_vz = vz[tid];
      rbmd::Real cont3 = 0.16666667;

      sum_vx +=  cont3 * (4 * fx[tid] - prev_fx[tid]) / mass[typei] * dt * fmt2v;
      sum_vy +=  cont3 * (4 * fy[tid] - prev_fy[tid]) / mass[typei] * dt * fmt2v;
      sum_vz +=  cont3 * (4 * fz[tid] - prev_fz[tid]) / mass[typei] * dt * fmt2v;

      sum_px += sum_vx * dt ;
      sum_py += sum_vy * dt ;
      sum_pz += sum_vz * dt ;

      // if (tid == 1) {
      //     printf("fx[2] = %f, fx_prev[2] = %f\n", fx[1], prev_fx[1]);
      // }
      prev_fx[tid] = fx[tid];
      prev_fy[tid] = fy[tid];
      prev_fz[tid] = fz[tid];
      // if (tid == 1) {
      //   printf("fx_prev[2] = %f\n", prev_fx[1]);
      // }

      px[tid] = sum_px;
      py[tid] = sum_py;
      pz[tid] = sum_pz;

      vx[tid] = sum_vx;
      vy[tid] = sum_vy;
      vz[tid] = sum_vz;

      ApplyPBC(box, px[tid], py[tid], pz[tid],
        flag_px[tid], flag_py[tid],flag_pz[tid]);
    }
  }
}

void UpdatePositionFlagOpbm<device::DEVICE_GPU>::operator()(
 const rbmd::Id num_atoms,const rbmd::Real fmt2v, const rbmd::Real dt, rbmd::Id test_current_step,Box  box  , const rbmd::Id* atoms_type,
 const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz,const rbmd::Real* mass,
 rbmd::Real* prev_fx, rbmd::Real* prev_fy, rbmd::Real* prev_fz,
 rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz,
 rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz,
 rbmd::Id* flag_px,rbmd::Id* flag_py, rbmd::Id* flag_pz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdatePositionFlagbm<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms,fmt2v, dt, test_current_step, box, atoms_type, fx, fy, fz, mass, prev_fx, prev_fy, prev_fz, vx, vy, vz, px, py, pz, flag_px, flag_py, flag_pz));
}

}  // namespace op
