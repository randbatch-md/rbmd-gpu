#include "update_temperature_op.h"
#include "rbmd_define.h"
#include <hipcub/hipcub.hpp>

namespace op
{
    #define THREADS_PER_BLOCK 256

	__global__
		void ComputeTemperature(const rbmd::Id num_atoms,
			                    const rbmd::Real mvv2e,
			                    const rbmd::Real* mass,
			                    const rbmd::Real* vx,
			                    const rbmd::Real* vy,
			                    const rbmd::Real* vz,
			                    rbmd::Real* temp_contrib)
	{
		__shared__ typename hipcub::BlockReduce<rbmd::Real, BLOCK_SIZE>::TempStorage temp_storage;

		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		rbmd::Real local_temp = 0;
		if (tid < num_atoms)
		{
			local_temp = mvv2e * mass[tid] * (vx[tid] * vx[tid] + vy[tid] * vy[tid] + vz[tid] * vz[tid]);
		}

		rbmd::Real block_sum = hipcub::BlockReduce<rbmd::Real, BLOCK_SIZE>(temp_storage).Sum(local_temp);
		if (threadIdx.x == 0) {
			atomicAdd(temp_contrib, block_sum);
		}

	}

	__global__
		void UpdataVelocityRescale(const rbmd::Id num_atoms,
			                       const rbmd::Real coeff_rescale,
			                       rbmd::Real* vx,
			                       rbmd::Real* vy,
			                       rbmd::Real* vz)
	{
		int tid = threadIdx.x + blockIdx.x * blockDim.x;

		if (tid < num_atoms)
		{
			vx[tid] = vx[tid] * coeff_rescale;
			vy[tid] = vy[tid] * coeff_rescale;
			vz[tid] = vz[tid] * coeff_rescale;
		}
	}

	__global__
		void UpdataVelocityNoseHoover(const rbmd::Id num_atoms,
			                          const rbmd::Real dt,
			                          const rbmd::Real fmt2v,
			                          const rbmd::Real nosehooverxi,
			                          const rbmd::Real* mass,
			                          const rbmd::Real* fx,
			                          const rbmd::Real* fy,
			                          const rbmd::Real* fz,
			                          rbmd::Real* vx,
			                          rbmd::Real* vy,
			                          rbmd::Real* vz)
	{
		int tid = threadIdx.x + blockIdx.x * blockDim.x;

		if (tid < num_atoms)
		{
			vx[tid] += 0.5 * dt * (fx[tid] / mass[tid] - nosehooverxi * vx[tid]) * fmt2v;
			vy[tid] += 0.5 * dt * (fy[tid] / mass[tid] - nosehooverxi * vy[tid]) * fmt2v;
			vz[tid] += 0.5 * dt * (fz[tid] / mass[tid] - nosehooverxi * vz[tid]) * fmt2v;
		}
	}

	__global__
		void UpdataVelocityBerendsen(const rbmd::Id num_atoms,
			                         const rbmd::Real coeff_Berendsen,
			                         rbmd::Real* vx,
			                         rbmd::Real* vy,
			                         rbmd::Real* vz)
	{
		int tid = threadIdx.x + blockIdx.x * blockDim.x;

		if (tid < num_atoms)
		{
			vx[tid] = vx[tid] * coeff_Berendsen;
			vy[tid] = vy[tid] * coeff_Berendsen;
			vz[tid] = vz[tid] * coeff_Berendsen;
		}
	}

	void ComputeTemperatureOp<device::DEVICE_GPU>::operator()(const rbmd::Id num_atoms,
			                                                       const rbmd::Real mvv2e,
			                                                       const rbmd::Real* mass,
			                                                       const rbmd::Real* vx,
			                                                       const rbmd::Real* vy,
			                                                       const rbmd::Real* vz,
			                                                       rbmd::Real* temp_contrib)
	{
		unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
		CHECK_KERNEL(ComputeTemperature <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>> (num_atoms, mvv2e, mass, vx, vy, vz, temp_contrib));
	}
	


	void UpdataVelocityRescaleOp<device::DEVICE_GPU>::operator()(const rbmd::Id num_atoms,
																	  const rbmd::Real coeff_rescale,
																	  rbmd::Real* vx,
																	  rbmd::Real* vy,
																	  rbmd::Real* vz)
	{
		unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
		CHECK_KERNEL(UpdataVelocityRescale <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>> (num_atoms, coeff_rescale, vx, vy, vz));
	}
	


	void UpdataVelocityNoseHooverOp<device::DEVICE_GPU>::operator()(const rbmd::Id num_atoms,
																		 const rbmd::Real dt,
																		 const rbmd::Real fmt2v,
																		 const rbmd::Real nosehooverxi,
																		 const rbmd::Real* mass,
																		 const rbmd::Real* fx,
																		 const rbmd::Real* fy,
																		 const rbmd::Real* fz,
																		 rbmd::Real* vx,
																		 rbmd::Real* vy,
																		 rbmd::Real* vz) 
	{
		unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
		CHECK_KERNEL(UpdataVelocityNoseHoover <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>> (num_atoms, dt, fmt2v, nosehooverxi, mass, fx, fy, fz, vx, vy, vz));
	}
	

	void UpdataVelocityBerendsenOp<device::DEVICE_GPU>::operator()(const rbmd::Id num_atoms,
			                                                            const rbmd::Real coeff_Berendsen,
			                                                            rbmd::Real* vx,
			                                                            rbmd::Real* vy,
			                                                            rbmd::Real* vz)
	{
		unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
		CHECK_KERNEL(UpdataVelocityBerendsen <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>> (num_atoms, coeff_Berendsen, vx, vy, vz));
	}
}

