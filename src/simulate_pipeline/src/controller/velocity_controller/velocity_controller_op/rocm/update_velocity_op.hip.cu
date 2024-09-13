#include "update_velocity_op.h"
#include "rbmd_define.h"

namespace op
{
    #define THREADS_PER_BLOCK 256

	__global__
		void UpdateVelocity(const rbmd::Id num_atoms,
			                const rbmd::Real dt,
			                const rbmd::Real fmt2v,
					        const rbmd::Id* atoms_type,
			                const rbmd::Real* fx,
			                const rbmd::Real* fy,
			                const rbmd::Real* fz,
			                const rbmd::Real* mass,
			                rbmd::Real* vx,
			                rbmd::Real* vy,
			                rbmd::Real* vz,
			                rbmd::Real* px,
			                rbmd::Real* py,
			                rbmd::Real* pz )
	{
		int tid = threadIdx.x + blockIdx.x * blockDim.x;

		if (tid < num_atoms)
		{
			rbmd::Id typei = atoms_type[tid]-1; 
			rbmd::Real massi = 1.0;

			//printf("--------test---massi:%f---\n",massi);

			// 更新速度
			//vx[tid] += 0.5 * fx[tid] / massi * dt * fmt2v;
			//vy[tid] += 0.5 * fy[tid] / massi * dt * fmt2v;
			//vz[tid] += 0.5 * fz[tid] / massi * dt * fmt2v;

			//// 更新位置
			//px[tid] += vx[tid] * dt;
			//py[tid] += vy[tid] * dt;
			//pz[tid] += vz[tid] * dt;


			rbmd::Real sum_vx = vx[tid];
			rbmd::Real sum_vy = vy[tid];
			rbmd::Real sum_vz = vz[tid];

			sum_vx += 0.5 * fx[tid] / massi * dt * fmt2v;
			sum_vy += 0.5 * fy[tid] / massi * dt * fmt2v;
			sum_vz += 0.5 * fz[tid] / massi * dt * fmt2v;

		    vx[tid] += sum_vx;
		    vy[tid] += sum_vy;
		    vz[tid] += sum_vz;

			//atomicAdd(&vx[tid], sum_vx);
			//atomicAdd(&vy[tid], sum_vy);
			//atomicAdd(&vz[tid], sum_vz);




			//
			rbmd::Real sum_px = px[tid];
			rbmd::Real sum_py = py[tid];
			rbmd::Real sum_pz = pz[tid];

			sum_px += sum_vx * dt;
			sum_py += sum_vy * dt;
			sum_pz += sum_vz * dt;

			px[tid] += sum_px;
			py[tid] += sum_py;
			pz[tid] += sum_pz;

			//atomicAdd(&px[tid], sum_px);
			//atomicAdd(&py[tid], sum_py);
			//atomicAdd(&pz[tid], sum_pz);


			//printf("--------test--vx:%f-,vy:%f,vz:%f--\n", vx[tid], vy[tid], vz[tid]);
			printf("--------test--px:%f-,py:%f,pz:%f--\n", px[tid], py[tid], pz[tid]);
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
		//printf("--------test---mass:%f---\n",mass[0]);
		unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
		CHECK_KERNEL(UpdateVelocity <<<blocks_per_grid, BLOCK_SIZE, 0, 0 >>> (num_atoms, dt, fmt2v, atoms_type,mass, fx, fy, fz, vx, vy, vz, px, py, pz));
	}
}

