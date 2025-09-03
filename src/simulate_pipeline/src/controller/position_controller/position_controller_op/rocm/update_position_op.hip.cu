#include "rbmd_define.h"
#include "update_position_op.h"

namespace op {
#define THREADS_PER_BLOCK 256

//蛙跳leapfrog
__global__ void UpdatePositionFlag(
const rbmd::Id num_atoms, const rbmd::Real dt, const rbmd::Real fmt2v, const rbmd::Id* atoms_type, const rbmd::Real* mass, Box  box,
const rbmd::Real* vx, const rbmd::Real* vy, const rbmd::Real* vz,
rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, rbmd::Id* flag_px,
rbmd::Id* flag_py, rbmd::Id* flag_pz, const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz) {
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

//leapfrog
void UpdatePositionFlagOp<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms, const rbmd::Real dt, const rbmd::Real fmt2v, const rbmd::Id* atoms_type, const rbmd::Real* mass, Box  box,
    const rbmd::Real* vx, const rbmd::Real* vy, const rbmd::Real* vz,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, rbmd::Id* flag_px,
    rbmd::Id* flag_py, rbmd::Id* flag_pz,const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdatePositionFlag<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms, dt, fmt2v, atoms_type, mass, box, vx, vy, vz, px, py, pz, flag_px, flag_py, flag_pz, fx, fy, fz));
}

// vv
 __global__ void UpdatePositionFlagvv(
 const rbmd::Id num_atoms, const rbmd::Real dt, const rbmd::Real fmt2v, const rbmd::Id* atoms_type, const rbmd::Real* mass, Box  box,
 const rbmd::Real* vx, const rbmd::Real* vy, const rbmd::Real* vz,
 rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, rbmd::Id* flag_px,
 rbmd::Id* flag_py, rbmd::Id* flag_pz, const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz) {
   int tid = threadIdx.x + blockIdx.x * blockDim.x;
   if (tid < num_atoms) {
     rbmd::Id typei = atoms_type[tid];
     rbmd::Real sum_px = px[tid];
     rbmd::Real sum_py = py[tid];
     rbmd::Real sum_pz = pz[tid];
     rbmd::Real const_half = 0.5000;

     sum_px += vx[tid] * dt + const_half * fx[tid] / mass[typei] * dt * dt * fmt2v;
     sum_py += vy[tid] * dt + const_half * fy[tid] / mass[typei] * dt * dt * fmt2v;
     sum_pz += vz[tid] * dt + const_half * fz[tid] / mass[typei] * dt * dt * fmt2v;

     // sum_px += vx[tid] * dt;
     // sum_py += vy[tid] * dt;
     // sum_pz += vz[tid] * dt;

     px[tid] = sum_px;
     py[tid] = sum_py;
     pz[tid] = sum_pz;

     ApplyPBC(box, px[tid], py[tid], pz[tid],
       flag_px[tid], flag_py[tid],flag_pz[tid]);
   }
}

// vv
void UpdatePositionFlagOpvv<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms, const rbmd::Real dt, const rbmd::Real fmt2v, const rbmd::Id* atoms_type, const rbmd::Real* mass, Box  box  ,
    const rbmd::Real* vx, const rbmd::Real* vy, const rbmd::Real* vz,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, rbmd::Id* flag_px,
    rbmd::Id* flag_py, rbmd::Id* flag_pz,const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdatePositionFlagvv<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms, dt, fmt2v, atoms_type, mass, box, vx, vy, vz, px, py, pz, flag_px, flag_py, flag_pz, fx, fy, fz));
}

// 4阶PRK
__global__ void UpdatePositionFlag1(
    const rbmd::Id num_atoms, const rbmd::Real d1, const rbmd::Real dt, const rbmd::Real fmt2v, const rbmd::Id* atoms_type, const rbmd::Real* mass, Box  box,
    const rbmd::Real* vx, const rbmd::Real* vy, const rbmd::Real* vz,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, rbmd::Id* flag_px,
    rbmd::Id* flag_py, rbmd::Id* flag_pz, const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    rbmd::Id typei = atoms_type[tid];
    rbmd::Real sum_px = px[tid];
    rbmd::Real sum_py = py[tid];
    rbmd::Real sum_pz = pz[tid];

    sum_px += d1 * vx[tid] * dt ;
    sum_py += d1 * vy[tid] * dt ;
    sum_pz += d1 * vz[tid] * dt ;

    px[tid] = sum_px;
    py[tid] = sum_py;
    pz[tid] = sum_pz;

    ApplyPBC(box, px[tid], py[tid], pz[tid],
      flag_px[tid], flag_py[tid],flag_pz[tid]);
  }
}
__global__ void UpdatePositionFlag2(
    const rbmd::Id num_atoms, const rbmd::Real d2, const rbmd::Real dt, const rbmd::Real fmt2v, const rbmd::Id* atoms_type, const rbmd::Real* mass, Box  box,
    const rbmd::Real* vx, const rbmd::Real* vy, const rbmd::Real* vz,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, rbmd::Id* flag_px,
    rbmd::Id* flag_py, rbmd::Id* flag_pz, const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    rbmd::Id typei = atoms_type[tid];
    rbmd::Real sum_px = px[tid];
    rbmd::Real sum_py = py[tid];
    rbmd::Real sum_pz = pz[tid];

    sum_px += d2 * vx[tid] * dt ;
    sum_py += d2 * vy[tid] * dt ;
    sum_pz += d2 * vz[tid] * dt ;

    px[tid] = sum_px;
    py[tid] = sum_py;
    pz[tid] = sum_pz;

    ApplyPBC(box, px[tid], py[tid], pz[tid],
      flag_px[tid], flag_py[tid],flag_pz[tid]);
  }
}
__global__ void UpdatePositionFlag3(
    const rbmd::Id num_atoms, const rbmd::Real d3, const rbmd::Real dt, const rbmd::Real fmt2v, const rbmd::Id* atoms_type, const rbmd::Real* mass, Box  box,
    const rbmd::Real* vx, const rbmd::Real* vy, const rbmd::Real* vz,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, rbmd::Id* flag_px,
    rbmd::Id* flag_py, rbmd::Id* flag_pz, const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    rbmd::Id typei = atoms_type[tid];
    rbmd::Real sum_px = px[tid];
    rbmd::Real sum_py = py[tid];
    rbmd::Real sum_pz = pz[tid];

    sum_px += d3 * vx[tid] * dt;
    sum_py += d3 * vy[tid] * dt;
    sum_pz += d3 * vz[tid] * dt;

    px[tid] = sum_px;
    py[tid] = sum_py;
    pz[tid] = sum_pz;

    ApplyPBC(box, px[tid], py[tid], pz[tid], flag_px[tid], flag_py[tid],
             flag_pz[tid]);
  }
}
__global__ void UpdatePositionFlag4(
    const rbmd::Id num_atoms, const rbmd::Real d4, const rbmd::Real dt, const rbmd::Real fmt2v, const rbmd::Id* atoms_type, const rbmd::Real* mass, Box  box,
    const rbmd::Real* vx, const rbmd::Real* vy, const rbmd::Real* vz,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, rbmd::Id* flag_px,
    rbmd::Id* flag_py, rbmd::Id* flag_pz, const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    rbmd::Id typei = atoms_type[tid];
    rbmd::Real sum_px = px[tid];
    rbmd::Real sum_py = py[tid];
    rbmd::Real sum_pz = pz[tid];

    sum_px += d4 * vx[tid] * dt ;
    sum_py += d4 * vy[tid] * dt ;
    sum_pz += d4 * vz[tid] * dt ;

    px[tid] = sum_px;
    py[tid] = sum_py;
    pz[tid] = sum_pz;

    ApplyPBC(box, px[tid], py[tid], pz[tid],
      flag_px[tid], flag_py[tid],flag_pz[tid]);
  }
}

//PRK
void UpdatePositionFlagOp1<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms, const rbmd::Real d1, const rbmd::Real dt, const rbmd::Real fmt2v, const rbmd::Id* atoms_type, const rbmd::Real* mass, Box  box  ,
    const rbmd::Real* vx, const rbmd::Real* vy, const rbmd::Real* vz,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, rbmd::Id* flag_px,
    rbmd::Id* flag_py, rbmd::Id* flag_pz,const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdatePositionFlag1<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms, d1, dt, fmt2v, atoms_type, mass, box, vx, vy, vz, px, py, pz, flag_px, flag_py, flag_pz, fx, fy, fz));
}
void UpdatePositionFlagOp2<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms, const rbmd::Real d2, const rbmd::Real dt, const rbmd::Real fmt2v, const rbmd::Id* atoms_type, const rbmd::Real* mass, Box  box  ,
    const rbmd::Real* vx, const rbmd::Real* vy, const rbmd::Real* vz,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, rbmd::Id* flag_px,
    rbmd::Id* flag_py, rbmd::Id* flag_pz,const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdatePositionFlag2<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms, d2, dt, fmt2v, atoms_type, mass, box, vx, vy, vz, px, py, pz, flag_px, flag_py, flag_pz, fx, fy, fz));
}
void UpdatePositionFlagOp3<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms, const rbmd::Real d3, const rbmd::Real dt, const rbmd::Real fmt2v, const rbmd::Id* atoms_type, const rbmd::Real* mass, Box  box  ,
    const rbmd::Real* vx, const rbmd::Real* vy, const rbmd::Real* vz,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, rbmd::Id* flag_px,
    rbmd::Id* flag_py, rbmd::Id* flag_pz,const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdatePositionFlag3<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms, d3, dt, fmt2v, atoms_type, mass, box, vx, vy, vz, px, py, pz, flag_px, flag_py, flag_pz, fx, fy, fz));
}
void UpdatePositionFlagOp4<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms, const rbmd::Real d4,const rbmd::Real dt, const rbmd::Real fmt2v, const rbmd::Id* atoms_type, const rbmd::Real* mass, Box  box  ,
    const rbmd::Real* vx, const rbmd::Real* vy, const rbmd::Real* vz,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, rbmd::Id* flag_px,
    rbmd::Id* flag_py, rbmd::Id* flag_pz,const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdatePositionFlag4<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms, d4, dt, fmt2v, atoms_type, mass, box, vx, vy, vz, px, py, pz, flag_px, flag_py, flag_pz, fx, fy, fz));
}

//未调用
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

//未调用
void UpdatePositionOp<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms, const rbmd::Real dt, Box  box  , const rbmd::Real* vx,
    const rbmd::Real* vy, const rbmd::Real* vz, rbmd::Real* px, rbmd::Real* py,
    rbmd::Real* pz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdatePosition<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms, dt, box, vx, vy, vz, px,py, pz));
}

// Beeman
__global__ void UpdatePositionFlagBeeman(
    const rbmd::Id num_atoms, const rbmd::Real dt,rbmd::Id test_current_step, const rbmd::Real fmt2v,
    const rbmd::Id* atoms_type, const rbmd::Real* mass, Box box,
    rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz,
    const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz,                 // F(t)
    rbmd::Real* f_pre1_x, rbmd::Real* f_pre1_y, rbmd::Real* f_pre1_z, // F(t-Δt)
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz,
    rbmd::Id* flag_px, rbmd::Id* flag_py, rbmd::Id* flag_pz)
{
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

      f_pre1_x[tid] = fx[tid];
      f_pre1_y[tid] = fy[tid];
      f_pre1_z[tid] = fz[tid];

      sum_vx +=  cont1 * f_pre1_x[tid] / mass[typei] * dt * fmt2v;
      sum_vy +=  cont1 * f_pre1_y[tid] / mass[typei] * dt * fmt2v;
      sum_vz +=  cont1 * f_pre1_z[tid] / mass[typei] * dt * fmt2v;

      px[tid] = sum_px;
      py[tid] = sum_py;
      pz[tid] = sum_pz;

      vx[tid] = sum_vx;
      vy[tid] = sum_vy;
      vz[tid] = sum_vz;

      ApplyPBC(box, px[tid], py[tid], pz[tid],
        flag_px[tid], flag_py[tid],flag_pz[tid]);
    }
    else{
      rbmd::Id typei = atoms_type[tid];
      rbmd::Real sum_px = px[tid];
      rbmd::Real sum_py = py[tid];
      rbmd::Real sum_pz = pz[tid];
      rbmd::Real sum_vx = vx[tid];
      rbmd::Real sum_vy = vy[tid];
      rbmd::Real sum_vz = vz[tid];
      rbmd::Real const_one_sixth = 0.16666667;


      sum_px += vx[tid] * dt + const_one_sixth * (4 * fx[tid] - f_pre1_x[tid]) / mass[typei] * dt * dt * fmt2v;
      sum_py += vy[tid] * dt + const_one_sixth * (4 * fy[tid] - f_pre1_y[tid]) / mass[typei] * dt * dt * fmt2v;
      sum_pz += vz[tid] * dt + const_one_sixth * (4 * fz[tid] - f_pre1_z[tid]) / mass[typei] * dt * dt * fmt2v;

      sum_vx +=  const_one_sixth * (5 * fx[tid] - f_pre1_x[tid]) / mass[typei] * dt * fmt2v;
      sum_vy +=  const_one_sixth * (5 * fy[tid] - f_pre1_y[tid]) / mass[typei] * dt * fmt2v;
      sum_vz +=  const_one_sixth * (5 * fz[tid] - f_pre1_z[tid]) / mass[typei] * dt * fmt2v;

      // sum_vx +=  const_one_sixth * (4 * fx[tid] - f_pre1_x[tid]) / mass[typei] * dt * fmt2v;
      // sum_vy +=  const_one_sixth * (4 * fy[tid] - f_pre1_y[tid]) / mass[typei] * dt * fmt2v;
      // sum_vz +=  const_one_sixth * (4 * fz[tid] - f_pre1_z[tid]) / mass[typei] * dt * fmt2v;

      f_pre1_x[tid] = fx[tid];
      f_pre1_y[tid] = fy[tid];
      f_pre1_z[tid] = fz[tid];

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
// ... 相应地创建 UpdatePositionOpBeeman ...
void UpdatePositionFlagOpBeeman<device::DEVICE_GPU>::operator()(
 const rbmd::Id num_atoms, const rbmd::Real dt, rbmd::Id test_current_step ,const rbmd::Real fmt2v,
 const rbmd::Id* atoms_type, const rbmd::Real* mass, Box box,
 rbmd::Real* vx, rbmd::Real* vy,rbmd::Real* vz,
 const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz,
 rbmd::Real* f_pre1_x,rbmd::Real* f_pre1_y,rbmd::Real* f_pre1_z,
 rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz,
 rbmd::Id* flag_px, rbmd::Id* flag_py, rbmd::Id* flag_pz)
{ unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdatePositionFlagBeeman<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
  num_atoms, dt,test_current_step, fmt2v, atoms_type, mass, box,
  vx, vy, vz,
  fx, fy, fz,
  f_pre1_x, f_pre1_y, f_pre1_z,
  px, py, pz, flag_px, flag_py, flag_pz))
}


}  // namespace op
