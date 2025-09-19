#include "rbmd_define.h"
#include "update_position_op.h"

namespace op {
#define THREADS_PER_BLOCK 256

__global__ void UpdatePositionFlag(
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

__global__ void UpdatePositionFlagbm(
    const rbmd::Id num_atoms,const rbmd::Real fmt2v, const rbmd::Real par_a,const rbmd::Real par_b, const rbmd::Real dt, rbmd::Id test_current_step,Box  box  , const rbmd::Id* atoms_type,
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
      rbmd::Real cont1 = 1.0/2.0;

      sum_px += vx[tid] * dt; // + cont1 * fx[tid] / mass[typei] * dt * dt * fmt2v;
      sum_py += vy[tid] * dt; // + cont1 * fy[tid] / mass[typei] * dt * dt * fmt2v;
      sum_pz += vz[tid] * dt; // + cont1 * fz[tid] / mass[typei] * dt * dt * fmt2v;

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
      rbmd::Real cont3 = 1.0/6.0;
      rbmd::Real cont4 = par_a/par_b;

      // sum_vx +=  cont3 * (4 * fx[tid] - prev_fx[tid]) / mass[typei] * dt * fmt2v;
      // sum_vy +=  cont3 * (4 * fy[tid] - prev_fy[tid]) / mass[typei] * dt * fmt2v;
      // sum_vz +=  cont3 * (4 * fz[tid] - prev_fz[tid]) / mass[typei] * dt * fmt2v;
      //
      // sum_px += sum_vx * dt ;
      // sum_py += sum_vy * dt ;
      // sum_pz += sum_vz * dt ;

      sum_px += sum_vx * dt + cont3 * (4 * fx[tid] - prev_fx[tid]) / mass[typei] * dt * dt * fmt2v;
      sum_py += sum_vy * dt + cont3 * (4 * fy[tid] - prev_fy[tid]) / mass[typei] * dt * dt * fmt2v;
      sum_pz += sum_vz * dt + cont3 * (4 * fz[tid] - prev_fz[tid]) / mass[typei] * dt * dt * fmt2v;

      sum_vx +=  ( (cont4 + cont3) * fx[tid] - cont3 * prev_fx[tid]) / mass[typei] * dt * fmt2v;
      sum_vy +=  ( (cont4 + cont3) * fy[tid] - cont3 * prev_fy[tid]) / mass[typei] * dt * fmt2v;
      sum_vz +=  ( (cont4 + cont3) * fz[tid] - cont3 * prev_fz[tid]) / mass[typei] * dt * fmt2v;


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

__global__ void UpdatePositionFlagbm1(
    const rbmd::Id num_atoms, const rbmd::Real fmt2v, const rbmd::Real par_a, const rbmd::Real par_b, const rbmd::Real dt, rbmd::Id test_current_step, Box box, const rbmd::Id* atoms_type,
    const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz, const rbmd::Real* mass,
     rbmd::Real* prev_fx,  rbmd::Real* prev_fy,  rbmd::Real* prev_fz, // Note: prev_fx/fy/fz are now const rbmd::Real*
    rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz,
    rbmd::Id* flag_px, rbmd::Id* flag_py, rbmd::Id* flag_pz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    rbmd::Id typei = atoms_type[tid];
    rbmd::Real inv_mass = 1.0 / mass[typei];
    rbmd::Real dt_sq = dt * dt * fmt2v;

    // Startup phase using Velocity-Verlet (first step)
    if (test_current_step < 3) {
      // This is the position update part of Velocity Verlet
      // r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt^2
      px[tid] += vx[tid] * dt  ;
      py[tid] += vy[tid] * dt ;
      pz[tid] += vz[tid] * dt ;

      // The first half of the velocity update for Velocity Verlet
      // v(t+dt/2) = v(t) + 0.5*a(t)*dt
      vx[tid] += 0.5 * fx[tid] * inv_mass * dt * fmt2v;
      vy[tid] += 0.5 * fy[tid] * inv_mass * dt * fmt2v;
      vz[tid] += 0.5 * fz[tid] * inv_mass * dt * fmt2v;

    }
    // Main loop using Beeman algorithm
    else {
      // Correct Beeman position update formula:
      // r(t+dt) = r(t) + v(t)*dt + (2/3)*a(t)*dt^2 - (1/6)*a(t-dt)*dt^2
      // Which is equivalent to: r(t) + v(t)*dt + (1/6)*(4*a(t) - a(t-dt))*dt^2
      const rbmd::Real cont3 = 0.1666667;
      px[tid] += vx[tid] * dt + cont3 * (4.0 * fx[tid] - prev_fx[tid]) * inv_mass * dt_sq;
      py[tid] += vy[tid] * dt + cont3 * (4.0 * fy[tid] - prev_fy[tid]) * inv_mass * dt_sq;
      pz[tid] += vz[tid] * dt + cont3 * (4.0 * fz[tid] - prev_fz[tid]) * inv_mass * dt_sq;


      // REMOVED: Incorrect velocity update was here.
      // REMOVED: Incorrect history update (prev_fx = fx) was here.
    }
    prev_fx[tid] = fx[tid];
    prev_fy[tid] = fy[tid];
    prev_fz[tid] = fz[tid];
    // Apply Periodic Boundary Conditions after position update
    ApplyPBC(box, px[tid], py[tid], pz[tid],
      flag_px[tid], flag_py[tid], flag_pz[tid]);
  }
}

__global__ void UpdatePositionFlagbm2(
    const rbmd::Id num_atoms, const rbmd::Real fmt2v, const rbmd::Real par_a, const rbmd::Real par_b, const rbmd::Real dt, rbmd::Id test_current_step, Box box, const rbmd::Id* atoms_type,
    const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz, const rbmd::Real* mass,
     rbmd::Real* prev_fx,  rbmd::Real* prev_fy,  rbmd::Real* prev_fz, // Note: prev_fx/fy/fz are now const rbmd::Real*
    rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz,
    rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz,
    rbmd::Id* flag_px, rbmd::Id* flag_py, rbmd::Id* flag_pz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    rbmd::Id typei = atoms_type[tid];
    rbmd::Real inv_mass = 1.0 / mass[typei];
    rbmd::Real dt_sq = dt * dt * fmt2v;

    // Startup phase using Velocity-Verlet (first step)
    if (test_current_step < 2) {
      // This is the position update part of Velocity Verlet
      // r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt^2
      px[tid] += vx[tid] * dt  ;
      py[tid] += vy[tid] * dt ;
      pz[tid] += vz[tid] * dt ;

      // The first half of the velocity update for Velocity Verlet
      // v(t+dt/2) = v(t) + 0.5*a(t)*dt
      vx[tid] += 0.5 * fx[tid] * inv_mass * dt * fmt2v;
      vy[tid] += 0.5 * fy[tid] * inv_mass * dt * fmt2v;
      vz[tid] += 0.5 * fz[tid] * inv_mass * dt * fmt2v;
    }
    // Main loop using Beeman algorithm
    else {
      // Which is equivalent to: r(t) + v(t)*dt + (1/2) * a(t) * t^2
      px[tid] += vx[tid] * dt +  0.5 * fx[tid] * inv_mass * dt_sq;
      py[tid] += vy[tid] * dt +  0.5 * fy[tid] * inv_mass * dt_sq;
      pz[tid] += vz[tid] * dt +  0.5 * fz[tid] * inv_mass * dt_sq;

      prev_fx[tid] = fx[tid];
      prev_fy[tid] = fy[tid];
      prev_fz[tid] = fz[tid];
    }

    // Apply Periodic Boundary Conditions after position update
    ApplyPBC(box, px[tid], py[tid], pz[tid],
      flag_px[tid], flag_py[tid], flag_pz[tid]);
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

void UpdatePositionFlagOpbm<device::DEVICE_GPU>::operator()(
 const rbmd::Id num_atoms,const rbmd::Real fmt2v, const rbmd::Real par_a,const rbmd::Real par_b, const rbmd::Real dt, rbmd::Id test_current_step,Box  box  , const rbmd::Id* atoms_type,
 const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz,const rbmd::Real* mass,
 rbmd::Real* prev_fx, rbmd::Real* prev_fy, rbmd::Real* prev_fz,
 rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz,
 rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz,
 rbmd::Id* flag_px,rbmd::Id* flag_py, rbmd::Id* flag_pz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdatePositionFlagbm<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms,fmt2v, par_a,par_b, dt, test_current_step, box, atoms_type, fx, fy, fz, mass, prev_fx, prev_fy, prev_fz, vx, vy, vz, px, py, pz, flag_px, flag_py, flag_pz));
}

}  // namespace op
