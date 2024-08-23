#include "temperature_controller_op/update_temperature_op.h"
#include "base/rocm.h"

namespace op
{
    #define THREADS_PER_BLOCK 256

	__global__
		void ComputeTemperature(const rbmd::Id& num_atoms,
			                    const rbmd::Real& mvv2e,
			                    const rbmd::Real* mass,
			                    const rbmd::Real* vx,
			                    const rbmd::Real* vy,
			                    const rbmd::Real* vz,
			                    rbmd::Real& temp_sum)
	{
		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		
		if (tid < num_atoms)
		{
			temp_sum += mvv2e * mass[tid] * (vx[tid] * vx[tid] + vy[tid] * vy[tid] + vz[tid] * vz[tid]);
		}
	}

	__global__
		void UpdataVelocity(const rbmd::Id& num_atoms,
			                const rbmd::Real& coeff_rescale,
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


	struct ComputeTemperatureOp<device::DEVICE_GPU>
	{
		void operator()(const rbmd::Id& num_atoms,
			            const rbmd::Real& mvv2e,
			            const rbmd::Real* mass,
			            const rbmd::Real* vx,
			            const rbmd::Real* vy,
			            const rbmd::Real* vz,
			            rbmd::Real& temp_sum)
		{
			int block_per_grid = (nAtoms + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

			hipLaunchKernelGGL(HIP_KERNEL_NAME(ComputeTemperature), dim3(block_per_grid), dim3(THREADS_PER_BLOCK), 0, 0,
				num_atoms, mvv2e, mass, vx, vy, vz, temp_sum);

			hipErrorCheck(hipGetLastError());
			hipErrorCheck(hipDeviceSynchronize());
		}
	};


	struct UpdataVelocityOp<device::DEVICE_GPU>
	{
		void operator()(const rbmd::Id& num_atoms,
			            const rbmd::Real& coeff_rescale,
			            rbmd::Real* vx,
			            rbmd::Real* vy,
			            rbmd::Real* vz)
		{
			int block_per_grid = (nAtoms + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

			hipLaunchKernelGGL(HIP_KERNEL_NAME(UpdataVelocity), dim3(block_per_grid), dim3(THREADS_PER_BLOCK), 0, 0,
				num_atoms, coeff_rescale, vx, vy, vz);

			hipErrorCheck(hipGetLastError());
			hipErrorCheck(hipDeviceSynchronize());
		}
	};


	template struct ComputeTemperatureOp<device::DEVICE_GPU>;
	template struct UpdataVelocityOp<device::DEVICE_GPU>;

}

