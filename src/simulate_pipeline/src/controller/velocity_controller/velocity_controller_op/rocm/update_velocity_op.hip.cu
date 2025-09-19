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

__global__ void UpdateVelocitybm(const rbmd::Id num_atoms, const rbmd::Real par_a, const rbmd::Real par_b,const rbmd::Real dt,
  rbmd::Id test_current_step,const rbmd::Real fmt2v,const rbmd::Id* atoms_type,
  const rbmd::Real* mass, const rbmd::Real* fx,const rbmd::Real* fy, const rbmd::Real* fz,
  rbmd::Real* prev_fx, rbmd::Real* prev_fy, rbmd::Real* prev_fz,
  rbmd::Real* pr_prev_fx, rbmd::Real* pr_prev_fy, rbmd::Real* pr_prev_fz,
  rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {

    if (test_current_step  < 2){
      rbmd::Id typei = atoms_type[tid];
      rbmd::Real sum_vx = vx[tid];
      rbmd::Real sum_vy = vy[tid];
      rbmd::Real sum_vz = vz[tid];
      rbmd::Real cont = 1.0/2.0;

      sum_vx +=  cont * fx[tid] / mass[typei] * dt * fmt2v;
      sum_vy +=  cont * fy[tid] / mass[typei] * dt * fmt2v;
      sum_vz +=  cont * fz[tid] / mass[typei] * dt * fmt2v;

      vx[tid] = sum_vx;
      vy[tid] = sum_vy;
      vz[tid] = sum_vz;

      pr_prev_fx[tid] = prev_fx[tid];
      pr_prev_fy[tid] = prev_fy[tid];
      pr_prev_fz[tid] = prev_fz[tid];

      prev_fx[tid] = fx[tid];
      prev_fy[tid] = fy[tid];
      prev_fz[tid] = fz[tid];
    }
    else {
      rbmd::Id typei = atoms_type[tid];
      rbmd::Real sum_vx = vx[tid];
      rbmd::Real sum_vy = vy[tid];
      rbmd::Real sum_vz = vz[tid];
      rbmd::Real cont3 = 1.0/3.0;
      rbmd::Real cont4 = par_a/par_b;


      sum_vx += (cont3 * fx[tid] + (cont3 * 2 - cont4) * prev_fx[tid]) / mass[typei] * dt * fmt2v;
      sum_vy += (cont3 * fy[tid] + (cont3 * 2 - cont4) * prev_fy[tid]) / mass[typei] * dt * fmt2v;
      sum_vz += (cont3 * fz[tid] + (cont3 * 2 - cont4) * prev_fz[tid]) / mass[typei] * dt * fmt2v;

      // sum_vx += (0.3333334 * fx[tid] + 0.2500000 * prev_fx[tid]) / mass[typei] * dt * fmt2v;
      // sum_vy += (0.3333334 * fy[tid] + 0.2500000 * prev_fy[tid]) / mass[typei] * dt * fmt2v;
      // sum_vz += (0.3333334 * fz[tid] + 0.2500000 * prev_fz[tid]) / mass[typei] * dt * fmt2v;

      // sum_vx += cont3 * (2 * fx[tid] + 5 * prev_fx[tid] - pr_prev_fx[tid]) / mass[typei] * dt * fmt2v;
      // sum_vy += cont3 * (2 * fy[tid] + 5 * prev_fy[tid] - pr_prev_fy[tid]) / mass[typei] * dt * fmt2v;
      // sum_vz += cont3 * (2 * fz[tid] + 5 * prev_fz[tid] - pr_prev_fz[tid]) / mass[typei] * dt * fmt2v;

      vx[tid] = sum_vx;
      vy[tid] = sum_vy;
      vz[tid] = sum_vz;

      pr_prev_fx[tid] = prev_fx[tid];
      pr_prev_fy[tid] = prev_fy[tid];
      pr_prev_fz[tid] = prev_fz[tid];

      prev_fx[tid] = fx[tid];
      prev_fy[tid] = fy[tid];
      prev_fz[tid] = fz[tid];
    }
  }
}

__global__ void UpdateVelocitybm1(const rbmd::Id num_atoms, const rbmd::Real par_a, const rbmd::Real par_b, const rbmd::Real dt,
  rbmd::Id test_current_step, const rbmd::Real fmt2v, const rbmd::Id* atoms_type,
  const rbmd::Real* mass, const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz,
  rbmd::Real* prev_fx, rbmd::Real* prev_fy, rbmd::Real* prev_fz,
  rbmd::Real* pr_prev_fx, rbmd::Real* pr_prev_fy, rbmd::Real* pr_prev_fz,
  rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    rbmd::Id typei = atoms_type[tid];
    rbmd::Real inv_mass = 1.0 / mass[typei];
    rbmd::Real dt_fmt = dt * fmt2v;

    // Startup phase using Velocity-Verlet
    if (test_current_step < 3) {
      // This is the second half of the velocity update for Velocity Verlet
      // v(t+dt) = v(t+dt/2) + 0.5*a(t+dt)*dt
      // Note: vx[tid] already holds v(t+dt/2) from the position kernel
      vx[tid] += 0.5 * fx[tid] * inv_mass * dt_fmt;
      vy[tid] += 0.5 * fy[tid] * inv_mass * dt_fmt;
      vz[tid] += 0.5 * fz[tid] * inv_mass * dt_fmt;
    }
    // Main loop using Beeman algorithm
    else {
      // CORRECTED: Implemented the canonical Beeman velocity formula.
      // v(t+dt) = v(t) + [ (1/3)a(t+dt) + (5/6)a(t) - (1/6)a(t-dt) ]*dt
      const rbmd::Real c1 = 0.33333333;
      const rbmd::Real c2 = 0.833333333;
      const rbmd::Real c3 = 0.1666666667;

      // Here: fx is a(t+dt), prev_fx is a(t), pr_prev_fx is a(t-dt)
      vx[tid] += (c1 * fx[tid] + c2 * prev_fx[tid] - c3 * pr_prev_fx[tid]) * inv_mass * dt_fmt;
      vy[tid] += (c1 * fy[tid] + c2 * prev_fy[tid] - c3 * pr_prev_fy[tid]) * inv_mass * dt_fmt;
      vz[tid] += (c1 * fz[tid] + c2 * prev_fz[tid] - c3 * pr_prev_fz[tid]) * inv_mass * dt_fmt;
    }

    // Update force history for the NEXT step AFTER velocity has been calculated
    // a(t-dt) for next step becomes a(t) from this step
    pr_prev_fx[tid] = prev_fx[tid];
    pr_prev_fy[tid] = prev_fy[tid];
    pr_prev_fz[tid] = prev_fz[tid];

    // a(t) for next step becomes a(t+dt) from this step
    prev_fx[tid] = fx[tid];
    prev_fy[tid] = fy[tid];
    prev_fz[tid] = fz[tid];
  }
}

__global__ void UpdateVelocitybm2(const rbmd::Id num_atoms, const rbmd::Real par_a, const rbmd::Real par_b, const rbmd::Real dt,
  rbmd::Id test_current_step, const rbmd::Real fmt2v, const rbmd::Id* atoms_type,
  const rbmd::Real* mass, const rbmd::Real* fx, const rbmd::Real* fy, const rbmd::Real* fz,
  rbmd::Real* prev_fx, rbmd::Real* prev_fy, rbmd::Real* prev_fz,
  rbmd::Real* pr_prev_fx, rbmd::Real* pr_prev_fy, rbmd::Real* pr_prev_fz,
  rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    rbmd::Id typei = atoms_type[tid];
    rbmd::Real inv_mass = 1.0 / mass[typei];
    rbmd::Real dt_fmt = dt * fmt2v;

    // Startup phase using Velocity-Verlet
    if (test_current_step < 2) {
      // This is the second half of the velocity update for Velocity Verlet
      // v(t+dt) = v(t+dt/2) + 0.5*a(t+dt)*dt
      // Note: vx[tid] already holds v(t+dt/2) from the position kernel
      vx[tid] += 0.5 * fx[tid] * inv_mass * dt_fmt;
      vy[tid] += 0.5 * fy[tid] * inv_mass * dt_fmt;
      vz[tid] += 0.5 * fz[tid] * inv_mass * dt_fmt;

    }
    else {
      //
      vx[tid] += 0.5 * (fx[tid] +  prev_fx[tid]) * inv_mass * dt_fmt;
      vy[tid] += 0.5 * (fy[tid] +  prev_fy[tid]) * inv_mass * dt_fmt;
      vz[tid] += 0.5 * (fz[tid] +  prev_fz[tid]) * inv_mass * dt_fmt;
    }

    // Update force history for the NEXT step AFTER velocity has been calculated
    // a(t-dt) for next step becomes a(t) from this step

    // a(t) for next step becomes a(t+dt) from this step
    prev_fx[tid] = fx[tid];
    prev_fy[tid] = fy[tid];
    prev_fz[tid] = fz[tid];
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

void UpdateVelocityOpbm<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms, const rbmd::Real par_a,const rbmd::Real par_b, const rbmd::Real dt, rbmd::Id test_current_step,
    const rbmd::Real fmt2v,const rbmd::Id* atoms_type,const rbmd::Real* mass,
    const rbmd::Real* fx,const rbmd::Real* fy, const rbmd::Real* fz,
    rbmd::Real* prev_fx, rbmd::Real* prev_fy, rbmd::Real* prev_fz,
    rbmd::Real* pr_prev_fx, rbmd::Real* pr_prev_fy, rbmd::Real* pr_prev_fz,
    rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdateVelocitybm<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms, par_a,par_b, dt, test_current_step, fmt2v, atoms_type, mass, fx, fy, fz, prev_fx, prev_fy, prev_fz, pr_prev_fx, pr_prev_fy, pr_prev_fz, vx, vy, vz));
}

}  // namespace op
