#include "rbmd_define.h"
#include "update_temperature_op.h"

namespace op {
#define THREADS_PER_BLOCK 256

__global__ void ComputeTemperature(const rbmd::Id num_atoms,
                                   const rbmd::Real mvv2e,
                                   const rbmd::Id* atoms_type,
                                   const rbmd::Real* mass, const rbmd::Real* vx,
                                   const rbmd::Real* vy, const rbmd::Real* vz,
                                   rbmd::Real* temp_contrib) {
  __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage
      temp_storage;

  rbmd::Real local_temp = 0;
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid < num_atoms) {
    local_temp = mvv2e * mass[atoms_type[tid]] *
                 (vx[tid] * vx[tid] + vy[tid] * vy[tid] + vz[tid] * vz[tid]);
  }

  rbmd::Real block_sum =
      BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>(temp_storage).Sum(local_temp);
  if (threadIdx.x == 0) {
    atomicAdd(temp_contrib, block_sum);
  }

}

__global__ void UpdataVelocityRescale(const rbmd::Id num_atoms,
                                      const rbmd::Real coeff_rescale,
                                      rbmd::Real* vx, rbmd::Real* vy,
                                      rbmd::Real* vz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    vx[tid] *= coeff_rescale;
    vy[tid] *= coeff_rescale;
    vz[tid] *= coeff_rescale;
  }
}

__global__ void UpdataVelocityNoseHoover(
    const rbmd::Id num_atoms, const rbmd::Real dt, const rbmd::Real fmt2v,
    const rbmd::Real nosehooverxi, const rbmd::Id* atoms_type,
    const rbmd::Real* mass, const rbmd::Real* fx, const rbmd::Real* fy,
    const rbmd::Real* fz, rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    vx[tid] += 0.5 * dt *
               (fx[tid] / mass[atoms_type[tid]] - nosehooverxi * vx[tid]) *
               fmt2v;
    vy[tid] += 0.5 * dt *
               (fy[tid] / mass[atoms_type[tid]] - nosehooverxi * vy[tid]) *
               fmt2v;
    vz[tid] += 0.5 * dt *
               (fz[tid] / mass[atoms_type[tid]] - nosehooverxi * vz[tid]) *
               fmt2v;
  }
}

__global__ void UpdataVelocityBerendsen(const rbmd::Id num_atoms,
                                        const rbmd::Real coeff_Berendsen,
                                        rbmd::Real* vx, rbmd::Real* vy,
                                        rbmd::Real* vz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    vx[tid] *= coeff_Berendsen;
    vy[tid] *= coeff_Berendsen;
    vz[tid] *= coeff_Berendsen;
  }
}

__global__
	void UpdataForceLangevin(const rbmd::Id num_atoms,
		                     const rbmd::Real gaussian_x,
		                     const rbmd::Real gaussian_y,
		                     const rbmd::Real gaussian_z,
		                     const rbmd::Real kbT,
		                     const rbmd::Real gamma,
		                     const rbmd::Real dt,
		                     const rbmd::Real* mass,
		                     const rbmd::Real* vx,
		                     const rbmd::Real* vy,
		                     const rbmd::Real* vz,
		                     rbmd::Real* fx,
		                     rbmd::Real* fy,
		                     rbmd::Real* fz)
	{
		int tid = threadIdx.x + blockIdx.x * blockDim.x;

		if (tid < num_atoms)
		{
			fx[tid] = fx[tid] - mass[tid] * gamma * vx[tid] + sqrt(2.0 * mass[tid] * gamma * kbT / dt) * gaussian_x;
			fy[tid] = fy[tid] - mass[tid] * gamma * vy[tid] + sqrt(2.0 * mass[tid] * gamma * kbT / dt) * gaussian_y;
			fz[tid] = fz[tid] - mass[tid] * gamma * vz[tid] + sqrt(2.0 * mass[tid] * gamma * kbT / dt) * gaussian_z;
		}
	}


// 新增的带COM移除的温度计算内核
__global__ void compute_temperature_com_kernel(const int num_atoms,
  const rbmd::Real mvv2e,const Real3 vbias,const int* atoms_type,
  const rbmd::Real* mass,const rbmd::Real* vx, const rbmd::Real* vy,
  const rbmd::Real* vz,rbmd::Real* temp_contrib)
{
  __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage
    temp_storage;

  rbmd::Real local_temp = 0;

  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid < num_atoms) {
    rbmd::Real dvx = vx[tid] - vbias.x;
    rbmd::Real dvy = vy[tid] - vbias.y;
    rbmd::Real dvz = vz[tid] - vbias.z;

    local_temp = mvv2e * mass[atoms_type[tid]] *
                 (dvx * dvx + dvy * dvy + dvz * dvz);
    //printf("local_temp: %f\n",local_temp);
  }

  rbmd::Real block_sum =
      BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>(temp_storage).Sum(local_temp);
  if (threadIdx.x == 0) {
    atomicAdd(temp_contrib, block_sum);
  }
}


// remove bias
__global__ void remove_bias_kernel(const int num_atoms,const Real3 vbias,
    rbmd:: Real* vx,rbmd::Real* vy,  rbmd::Real* vz)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < num_atoms) {
    // Remove bias,
    vx[i] = vx[i] - vbias.x;
    vy[i] = vy[i] - vbias.y;
    vz[i] = vz[i] - vbias.z;
  }
}

__global__ void restore_bias_kernel(const int num_atoms,const Real3 vbias,
    rbmd:: Real* vx,rbmd::Real* vy,  rbmd::Real* vz)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < num_atoms) {
    // restore bias,
    vx[i] = vx[i] + vbias.x;
    vy[i] = vy[i] + vbias.y;
    vz[i] = vz[i] + vbias.z;
  }
}
void ComputeTemperatureOp<device::DEVICE_GPU>::operator()(const rbmd::Id num_atoms,
		                                                  const rbmd::Real mvv2e,
								 const rbmd::Id* atoms_type,
		                                                  const rbmd::Real* mass,
		                                                  const rbmd::Real* vx,
		                                                  const rbmd::Real* vy,
		                                                  const rbmd::Real* vz,
		                                                  rbmd::Real* temp_contrib)

{
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(ComputeTemperature<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms, mvv2e, atoms_type, mass, vx, vy, vz, temp_contrib));
}

void UpdataVelocityRescaleOp<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms, const rbmd::Real coeff_rescale, rbmd::Real* vx,
    rbmd::Real* vy, rbmd::Real* vz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdataVelocityRescale<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms, coeff_rescale, vx, vy, vz));
}

void UpdataVelocityNoseHooverOp<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms, const rbmd::Real dt, const rbmd::Real fmt2v,
    const rbmd::Real nosehooverxi, const rbmd::Id* atoms_type,
    const rbmd::Real* mass, const rbmd::Real* fx, const rbmd::Real* fy,
    const rbmd::Real* fz, rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdataVelocityNoseHoover<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms, dt, fmt2v, nosehooverxi, atoms_type, mass, fx, fy, fz, vx, vy,
      vz));
}

void UpdataVelocityBerendsenOp<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms, const rbmd::Real coeff_Berendsen, rbmd::Real* vx,
    rbmd::Real* vy, rbmd::Real* vz) {
    unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
    CHECK_KERNEL(UpdataVelocityBerendsen <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>> (
        num_atoms, coeff_Berendsen, vx, vy, vz));
}

void UpdataForceLangevinOp<device::DEVICE_GPU>::operator()(const rbmd::Id num_atoms,
	                                                       const rbmd::Real gaussian_x,
	                                                       const rbmd::Real gaussian_y,
	                                                       const rbmd::Real gaussian_z,
	                                                       const rbmd::Real kbT,
	                                                       const rbmd::Real gamma,
	                                                       const rbmd::Real dt,
	                                                       const rbmd::Real* mass,
	                                                       const rbmd::Real* vx,
	                                                       const rbmd::Real* vy,
	                                                       const rbmd::Real* vz,
	                                                       rbmd::Real* fx,
	                                                       rbmd::Real* fy,
	                                                       rbmd::Real* fz)
	{
		unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
		CHECK_KERNEL(UpdataForceLangevin <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>> (num_atoms, gaussian_x, gaussian_y, gaussian_z, kbT, gamma, dt, mass, vx, vy, vz, fx, fy, fz));
	}

void ComputeTemperatureCOMOp<device::DEVICE_GPU>::operator()(const int num_atoms,
  const rbmd::Real mvv2e,const Real3 vbias,const int* atoms_type,
  const rbmd::Real* mass,const rbmd::Real* vx, const rbmd::Real* vy,
  const rbmd::Real* vz,rbmd::Real* temp_contrib) {
  // ... 计算blocks和threads ...
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(compute_temperature_com_kernel<<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>>
    (num_atoms,mvv2e,vbias,atoms_type,mass,vx,vy,vz,temp_contrib));
}


void RemoveBiasOp<device::DEVICE_GPU>::operator()(const int num_atoms,
  const Real3 vbias,rbmd:: Real* vx,rbmd::Real* vy,rbmd::Real* vz) {

  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(remove_bias_kernel<<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>>
    (num_atoms,vbias,vx,vy,vz));
}

void RestoreBiasOp<device::DEVICE_GPU>::operator()(const int num_atoms,
  const Real3 vbias,rbmd:: Real* vx,rbmd::Real* vy,  rbmd::Real* vz) {

  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(restore_bias_kernel<<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>>
    (num_atoms,vbias,vx,vy,vz));
}

}  // namespace op
