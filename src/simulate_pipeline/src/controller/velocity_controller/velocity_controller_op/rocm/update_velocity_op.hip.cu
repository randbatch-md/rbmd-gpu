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
			                rbmd::Real* vz,
			                rbmd::Real* px,
			                rbmd::Real* py,
			                rbmd::Real* pz )
	{
		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		//printf("--force_out----mass[0]---0:%f\n", mass[0]);
		//printf("--force_out----mass[0]---0:%d\n", num_atoms);

		if (tid < num_atoms)
		{
			rbmd::Id typei = atoms_type[tid] - 1;
			vx[tid] += 0.5 * fx[tid] / mass[typei] * dt * fmt2v;
			vy[tid] += 0.5 * fy[tid] / mass[typei] * dt * fmt2v;
			vz[tid] += 0.5 * fz[tid] / mass[typei] * dt * fmt2v;

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
			                                              rbmd::Real* vz,
														  rbmd::Real* px,
															rbmd::Real* py,
															rbmd::Real* pz)
	{
<<<<<<< .mine
		//printf("--------test---mass:%f---\n",mass[0]);


















=======
		//printf("---force_out----num_atoms:%d\n", num_atoms);
		//printf("---force_out----dt:%f\n", dt);
		//printf("---force_out----fmt2v:%f\n", fmt2v);
		//printf("---force_out----atoms_type:%f\n", atoms_type[0]);
		//printf("---force_out----mass[0]---0:%f\n", mass[0]);
		//printf("---force_out----mass[1]---0:%f\n", mass[0]);

		//printf("---force_out----fx:%f\n", fx[0]);
		//printf("---force_out----fy:%f\n", fy[0]);
		//printf("---force_out----fz:%f\n", fz[0]);
		//printf("---force_out----vx:%f\n", vx[0]);
		//printf("---force_out----vy:%f\n", vy[0]);
		//printf("---force_out----vz:%f\n", vz[0]);

		//for (size_t i = 0; i < num_atoms; i++)
		//{
		//	printf("---force_out----fx:%f\n", fx[i]);
		//}

>>>>>>> .theirs
		unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
<<<<<<< .mine
		CHECK_KERNEL(UpdateVelocity <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>> (num_atoms, dt, fmt2v, atoms_type,mass, fx, fy, fz, vx, vy, vz, px, py, pz));




=======
 		CHECK_KERNEL(UpdateVelocity <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>> (num_atoms, dt, fmt2v, atoms_type,mass, fx, fy, fz, vx, vy, vz));
		//for (size_t i = 0; i < num_atoms; i++)
		//{
		//	printf("---force_out----v:[ %f,%f,%f ]\n", vx[i], vy[i], vz[i]);
		//}
>>>>>>> .theirs
	}
}

