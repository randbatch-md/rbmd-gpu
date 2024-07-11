#include "near_force/direct_truncation/direct_truncation_op.h"
#include "base/rocm.h"

namespace op
{

#define THREADS_PER_BLOCK 256

template<typename FPTYPE>
__device__
void ComputeNeighbors(
	const rbmd::Real3* pos,
	const int& num_particles,
	const float& cut_off,
	int* neighborList)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x ;
	if (idx >= num_particles)
	{
		return;
	}

	auto p = pos[idx];
	int neighborCount = 0;

	for (int i = 0; i < num_particles; ++i)
	{
		if (i == idx)
		{
			continue;
		}

		auto q = pos[i];


	}



}


template <typename FPTYPE>
__device__ __inline__
void UpdateVelocity()
{

}

template<typename FPTYPE>
__global__ 
void ComputeForce(
	const FPTYPE* dt,
	const FPTYPE* fmt2v,
	const FPTYPE* mass,
	const Locator* locator,
	rbmd::Real3* position,
	rbmd::Real3* v,
	rbmd::Real3* force)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (100 == tid)
	{
		auto cell_id = locator->GetCellId(position[tid]);
		printf("cell id: %d %d %d\n", cell_id.x, cell_id.y, cell_id.z);
		printf("dt: %f\n", *dt);
		printf("fmt2v: %f\n", *fmt2v);
		printf("mass: %f\n", mass[0]);
		printf("v: %f\n", v[0].data[0]);
		printf("force: %f\n", force[0].data[0]);
	}
}


template<typename FPTYPE>
struct direct_truncation_op<FPTYPE, device::DEVICE_GPU>
{
	void operator()(
		const int& nSteps,
		const int& nAtoms,
		const rbmd::Real3& left,
		const rbmd::Real3& right,
		const rbmd::Id3& dim,
		const FPTYPE* dt,
		const FPTYPE* fmt2v,
		const FPTYPE* mass,
		const Locator* locator,
		rbmd::Real3* position,
		rbmd::Real3* v,
		rbmd::Real3* force)
	{
		int block = (nAtoms + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

		printf("nSteps: %d\n", nSteps);
		printf("nAtoms: %d\n", nAtoms);

		//cell id list
		int cell_ids[nAtoms];


		hipLaunchKernelGGL(HIP_KERNEL_NAME(ComputeForce<FPTYPE>), dim3(block), dim3(THREADS_PER_BLOCK), 0, 0,
			dt, fmt2v, mass, locator, position, v, force);


		hipErrorCheck(hipGetLastError());
		hipErrorCheck(hipDeviceSynchronize());
	}
};

template struct direct_truncation_op<float, device::DEVICE_GPU>;
template struct direct_truncation_op<double, device::DEVICE_GPU>;
}

