#include "update_velocity_op.h"
#include "rbmd_define.h"

namespace op
{
    #define THREADS_PER_BLOCK 256

	__global__
		void UpdateVelocity(const rbmd::Id& num_atoms,
			                const rbmd::Real& dt,
			                const rbmd::Real& fmt2v,
			                const rbmd::Real* fx,
			                const rbmd::Real* fy,
			                const rbmd::Real* fz,
			                const rbmd::Real* mass,
			                rbmd::Real* vx,
			                rbmd::Real* vy,
			                rbmd::Real* vz)
	{
		int tid = threadIdx.x + blockIdx.x * blockDim.x;

		if (tid < num_atoms)
		{
			vx[tid] += 0.5 * fx[tid] / mass[tid] * dt * fmt2v;
			vy[tid] += 0.5 * fy[tid] / mass[tid] * dt * fmt2v;
			vz[tid] += 0.5 * fz[tid] / mass[tid] * dt * fmt2v;
		}
	}


	void UpdateVelocityOp<device::DEVICE_GPU>::operator()(const rbmd::Id& num_atoms,
			                                              const rbmd::Real& dt,
			                                              const rbmd::Real& fmt2v,
			                                              const rbmd::Real* mass,
			                                              const rbmd::Real* fx,
			                                              const rbmd::Real* fy,
			                                              const rbmd::Real* fz,
			                                              rbmd::Real* vx,
			                                              rbmd::Real* vy,
			                                              rbmd::Real* vz)
	{
		unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
		CHECK_KERNEL(UpdateVelocity <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>> (num_atoms, dt, fmt2v, mass, fx, fy, fz, vx, vy, vz));
	}
}

