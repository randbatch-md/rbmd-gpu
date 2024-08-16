#include "velocity_controller/rocm/update_velocity_op.h"
#include "base/rocm.h"

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

		if (tid < nAtoms)
		{
			vx += 0.5 * fx / mass * dt * fmt2v;
			vy += 0.5 * fy / mass * dt * fmt2v;
			vz += 0.5 * fz / mass * dt * fmt2v;
		}
	}


	struct UpdateVelocityOp<device::DEVICE_GPU>
	{
		void operator()(const rbmd::Id& num_atoms,
			const rbmd::Real& dt,
			const rbmd::Real& fmt2v,
			const std::vector<rbmd::Real>* mass,
			const std::vector<rbmd::Real>* fx,
			const std::vector<rbmd::Real>* fy,
			const std::vector<rbmd::Real>* fz,
			std::vector<rbmd::Real>* vx,
			std::vector<rbmd::Real>* vy,
			std::vector<rbmd::Real>* vz);
		{
			int block = (nAtoms + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

			hipLaunchKernelGGL(HIP_KERNEL_NAME(UpdateVelocity), dim3(block), dim3(THREADS_PER_BLOCK), 0, 0,
				num_atoms, dt, fmt2v, mass, fx, fy, fz, vx, vy, vz);

			hipErrorCheck(hipGetLastError());
			hipErrorCheck(hipDeviceSynchronize());
		}
	};

	template struct UpdateVelocityOp<device::DEVICE_GPU>;
	template struct UpdateVelocityOp<device::DEVICE_GPU>;
}

