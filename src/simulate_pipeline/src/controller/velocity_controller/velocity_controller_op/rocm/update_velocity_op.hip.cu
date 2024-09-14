#include "update_velocity_op.h"
#include "rbmd_define.h"
#include <hip/hip_runtime.h>
namespace op
{
    #define THREADS_PER_BLOCK 256

	__global__
		void UpdateVelocity(const rbmd::Id num_atoms,
			                const rbmd::Real dt,
			                const rbmd::Real fmt2v,
					        const rbmd::Id* atoms_type,
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
			rbmd::Id typei = atoms_type[tid] - 1;
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


	void UpdateVelocityOp<device::DEVICE_GPU>::operator()(const rbmd::Id num_atoms,
			                                              const rbmd::Real dt,
			                                              const rbmd::Real fmt2v,
														  const rbmd::Id* atoms_type,
			                                              const rbmd::Real* mass,
			                                              const rbmd::Real* fx,
			                                              const rbmd::Real* fy,
			                                              const rbmd::Real* fz,
			                                              rbmd::Real* vx,
			                                              rbmd::Real* vy,
			                                              rbmd::Real* vz)
	{
		unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
 		CHECK_KERNEL(UpdateVelocity <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>> (num_atoms, dt, fmt2v, atoms_type,mass, fx, fy, fz, vx, vy, vz));

	}
}

