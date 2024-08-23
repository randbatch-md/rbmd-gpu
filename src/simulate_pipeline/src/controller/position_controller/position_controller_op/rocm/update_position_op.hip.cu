#include "position_controller_op/update_position_op.h"
#include "base/rocm.h"

namespace op
{
#define THREADS_PER_BLOCK 256

	__device__
		void UpdateFlagOverRangePoint(const rbmd::Real& min_x_tid,
			                          const rbmd::Real& min_y_tid,
			                          const rbmd::Real& min_z_tid,
			                          const rbmd::Real& max_x_tid,
			                          const rbmd::Real& max_y_tid,
			                          const rbmd::Real& max_z_tid,
                                      rbmd::Real& px_tid,
			                          rbmd::Real& py_tid,
			                          rbmd::Real& pz_tid,
			                          rbmd::Id& flag_px_tid,
			                          rbmd::Id& flag_py_tid,
			                          rbmd::Id& flag_pz_tid)
	{

		flag_px_tid += (px_tid > max_x_tid)-(px_tid < min_x_tid);
		px_tid += (px_tid < min_x_tid) * (max_x_tid - min_x_tid) - (px_tid > max_x_tid) * (max_x_tid - min_x_tid);

		flag_py_tid += (py_tid > max_y_tid) - (py_tid < min_y_tid);
		py_tid += (py_tid < min_y_tid) * (max_y_tid - min_y_tid) - (py_tid > max_y_tid) * (max_y_tid - min_y_tid);

		flag_pz_tid += (pz_tid > max_z_tid) - (pz_tid < min_z_tid);
		pz_tid += (pz_tid < min_z_tid) * (max_z_tid - min_z_tid) - (pz_tid > max_z_tid) * (max_z_tid - min_z_tid);
		/*
		if (px_tid < min_x_tid)
		{
			px_tid += max_x_tid - min_x_tid;
			flag_px_tid -= 1;
		}
		else if (px_tid > max_x_tid)
		{
			px_tid -= max_x_tid - min_x_tid;
			flag_px_tid += 1;
		}
		
		if (py_tid < min_y_tid)
		{
			py_tid += max_y_tid - min_y_tid;
			flag_py_tid -= 1;
		}
		else if (py_tid > max_y_tid)
		{
			py_tid -= max_y_tid - min_y_tid;
			flag_py_tid += 1;
		}

		if (pz_tid < min_z_tid)
		{
			pz_tid += max_z_tid - min_z_tid;
			flag_pz_tid -= 1;
		}
		else if (pz_tid > max_z_tid)
		{
			pz_tid -= max_z_tid - min_z_tid;
			flag_pz_tid += 1;
		}
		*/
	}

	__device__
		void UpdateOverRangePoint(const rbmd::Real& min_x_tid,
		                          const rbmd::Real& min_y_tid,
		                          const rbmd::Real& min_z_tid,
		                          const rbmd::Real& max_x_tid,
		                          const rbmd::Real& max_y_tid,
		                          const rbmd::Real& max_z_tid,
		                          rbmd::Real& px_tid,
		                          rbmd::Real& py_tid,
		                          rbmd::Real& pz_tid)
	{
		px_tid += (px_tid < min_x_tid) * (max_x_tid - min_x_tid) - (px_tid > max_x_tid) * (max_x_tid - min_x_tid);

		py_tid += (py_tid < min_y_tid) * (max_y_tid - min_y_tid) - (py_tid > max_y_tid) * (max_y_tid - min_y_tid);

		pz_tid += (pz_tid < min_z_tid) * (max_z_tid - min_z_tid) - (pz_tid > max_z_tid) * (max_z_tid - min_z_tid);
		
		/*
		if (px_tid < min_x_tid)
		{
			px_tid += max_x_tid - min_x_tid;
		}
		else if (px_tid > max_x_tid)
		{
			px_tid -= max_x_tid - min_x_tid;
		}

		if (py_tid < min_y_tid)
		{
			py_tid += max_y_tid - min_y_tid;
		}
		else if (py_tid > max_y_tid)
		{
			py_tid -= max_y_tid - min_y_tid;
		}

		if (pz_tid < min_z_tid)
		{
			pz_tid += max_z_tid - min_z_tid;
		}
		else if (pz_tid > max_z_tid)
		{
			pz_tid -= max_z_tid - min_z_tid;
		}
		*/
	}

	__global__
		void UpdatePositionFlag(const rbmd::Id& num_atoms,
			                    const rbmd::Real& dt,
			                    const rbmd::Real& min_x,
			                    const rbmd::Real& min_y,
			                    const rbmd::Real& min_z,
			                    const rbmd::Real& max_x,
			                    const rbmd::Real& max_y,
			                    const rbmd::Real& max_z,
			                    const rbmd::Real* vx,
			                    const rbmd::Real* vy,
			                    const rbmd::Real* vz,
			                    rbmd::Real* px,
			                    rbmd::Real* py,
			                    rbmd::Real* pz,
			                    rbmd::Id* flag_px,
			                    rbmd::Id* flag_py,
			                    rbmd::Id* flag_pz)
	{
		int tid = threadIdx.x + blockIdx.x * blockDim.x;

		if (tid < num_atoms)
		{
			px[tid] += vx[tid] * dt;
			py[tid] += vy[tid] * dt;
			pz[tid] += vz[tid] * dt;

			UpdateFlagOverRangePoint(min_x[tid],
				                     min_y[tid],
				                     min_z[tid],
				                     max_x[tid],
				                     max_y[tid],
				                     max_z[tid],
				                     px[tid],
				                     py[tid],
				                     pz[tid],
				                     flag_px[tid],
				                     flag_py[tid],
				                     flag_pz[tid]);
		}
	}

	__global__
		void UpdatePosition(const rbmd::Id& num_atoms,
			                const rbmd::Real& dt,
			                const rbmd::Real& min_x,
			                const rbmd::Real& min_y,
			                const rbmd::Real& min_z,
			                const rbmd::Real& max_x,
			                const rbmd::Real& max_y,
			                const rbmd::Real& max_z,
			                const rbmd::Real* vx,
			                const rbmd::Real* vy,
			                const rbmd::Real* vz,
			                rbmd::Real* px,
			                rbmd::Real* py,
			                rbmd::Real* pz)
	{
		int tid = threadIdx.x + blockIdx.x * blockDim.x;

		if (tid < num_atoms)
		{
			px[tid] += vx[tid] * dt;
			py[tid] += vy[tid] * dt;
			pz[tid] += vz[tid] * dt;

			UpdateOverRangePoint(min_x[tid],
				                 min_y[tid],
				                 min_z[tid],
				                 max_x[tid],
				                 max_y[tid],
				                 max_z[tid],
				                 px[tid],
				                 py[tid],
				                 pz[tid]);
		}
	}


	struct UpdatePositionFlagOp<device::DEVICE_GPU>
	{
		void operator()(const rbmd::Id& num_atoms,
			            const rbmd::Real& dt,
			            const rbmd::Real& min_x,
			            const rbmd::Real& min_y,
			            const rbmd::Real& min_z,
			            const rbmd::Real& max_x,
			            const rbmd::Real& max_y,
			            const rbmd::Real& max_z,
			            const rbmd::Real* vx,
			            const rbmd::Real* vy,
			            const rbmd::Real* vz,
			            rbmd::Real* px,
			            rbmd::Real* py,
			            rbmd::Real* pz,
			            rbmd::Id* flag_px,
			            rbmd::Id* flag_py,
			            rbmd::Id* flag_pz)
		{
			int block_per_grid = (nAtoms + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

			hipLaunchKernelGGL(HIP_KERNEL_NAME(UpdatePositionFlag), dim3(block_per_grid), dim3(THREADS_PER_BLOCK), 0, 0,
				num_atoms, dt, min_x, min_y, min_z, max_x, max_y, max_z, vx, vy, vz, px, py, pz, flag_px, flag_py, flag_pz);

			hipErrorCheck(hipGetLastError());
			hipErrorCheck(hipDeviceSynchronize());
		}
	};

	struct UpdatePositionOp<device::DEVICE_GPU>
	{
		void operator()(const rbmd::Id& num_atoms,
			            const rbmd::Real& dt,
			            const rbmd::Real& min_x,
			            const rbmd::Real& min_y,
			            const rbmd::Real& min_z,
			            const rbmd::Real& max_x,
			            const rbmd::Real& max_y,
			            const rbmd::Real& max_z,
			            const rbmd::Real* vx,
			            const rbmd::Real* vy,
			            const rbmd::Real* vz,
			            rbmd::Real* px,
			            rbmd::Real* py,
			            rbmd::Real* pz)
		{
			int block_per_grid = (nAtoms + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

			hipLaunchKernelGGL(HIP_KERNEL_NAME(UpdatePosition), dim3(block_per_grid), dim3(THREADS_PER_BLOCK), 0, 0,
				num_atoms, dt, min_x, min_y, min_z, max_x, max_y, max_z, vx, vy, vz, px, py, pz, flag_px, flag_py, flag_pz);

			hipErrorCheck(hipGetLastError());
			hipErrorCheck(hipDeviceSynchronize());
		}
	};

	template struct UpdatePositionFlagOp<device::DEVICE_GPU>;
	template struct UpdatePositionOp<device::DEVICE_GPU>;
}