#include "near_force/direct_truncation/direct_truncation_op.h"
#include "base/rocm.h"

namespace op
{

#define THREADS_PER_BLOCK 256

template <typename FPTYPE>
__device__ __inline__
void UpdateVelocity()
{

}

template<typename FPTYPE>
__global__ void ComputeForce(
	FPTYPE* dt,
	FPTYPE* fmt2v,
	FPTYPE* mass,
	FPTYPE* v,
	FPTYPE* force)
{
	printf("dt: %f\n", dt);
	printf("fmt2v: %f\n", fmt2v);
}

template<typename FPTYPE>
struct direct_truncation_op<FPTYPE, device::DEVICE_GPU>
{
	void operator()(
		const int& nSteps,
		const int& nAtoms,
		const FPTYPE* dt,
		const FPTYPE* fmt2v,
		const FPTYPE* mass,
		FPTYPE* v,
		FPTYPE* force)
	{
		int block = (nAtoms + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

		printf("nSteps: %d\n", nSteps);
		printf("nAtoms: %d\n", nAtoms);

		hipLaunchKernelGGL(HIP_KERNEL_NAME(ComputeForce<FPTYPE>), dim3(block), dim3(THREADS_PER_BLOCK), 0, 0,
			dt, fmt2v, mass, v, force);

		hipErrorCheck(hipGetLastError());
		hipErrorCheck(hipDeviceSynchronize());
	}
};

template struct direct_truncation_op<float, device::DEVICE_GPU>;
template struct direct_truncation_op<double, device::DEVICE_GPU>;
}

