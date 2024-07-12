#include "near_force/direct_truncation/direct_truncation_op.h"
#include "base/rocm.h"

namespace op
{

#define THREADS_PER_BLOCK 256

__device__
void ComputeCellId(
	const rbmd::Real3& position,
	rbmd::Id3& cellids,
	const rbmd::Real3& left,
	const rbmd::Real3& right,
	const rbmd::Id3& dim)
{
	printf("right: %f,%f,%f\n", right.data[0], right.data[1], right.data[2]);
	printf("left:%f,%f,%f\n", left.data[0], left.data[1], left.data[2]);
	printf("dim: %d,%d,%d\n", dim.data[0], dim.data[1], dim.data[2]);

	rbmd::Real3 dxdydz = (right - left) / dim; //should be shared memory
	printf("dxdydz: %f,%f,%f\n", dxdydz.data[0], dxdydz.data[1], dxdydz.data[2]);

	cellids.data[0] = (position.data[0] - left.data[0]) / dxdydz.data[0];
	cellids.data[1] = (position.data[1] - left.data[1]) / dxdydz.data[1];
	cellids.data[2] = (position.data[2] - left.data[2]) / dxdydz.data[2];
	printf("dim: %d,%d,%d", cellids.data[0], cellids.data[1], cellids.data[2]);

}

template <typename FPTYPE>
__device__ __inline__
void UpdateVelocity()
{

}

template<typename FPTYPE>
__global__ 
void ComputeForce(
	int nAtoms,
	FPTYPE* dt,
	FPTYPE* fmt2v,
	FPTYPE* mass,
	rbmd::Real3 left,
	rbmd::Real3 right,
	rbmd::Id3& dim,
	rbmd::Id3* cellid,
	//const Locator* locator,
	rbmd::Real3* position,
	rbmd::Real3* v,
	rbmd::Real3* force)
{
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid > nAtoms)
	{
		return;
	}

	//cell id list
	if (100 == tid)
	{
		ComputeCellId(position[tid], cellid[tid], left, right, dim);

		printf("cell id: %d,%d,%d", cellid[tid].data[0], cellid[tid].data[1], cellid[tid].data[2]);
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
		const rbmd::Real3& left, //shared memory
		const rbmd::Real3& right,
		const rbmd::Id3& dim,
		rbmd::Id3* cellid,
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
		printf("left: %f\n", left.data[0]);
		printf("right: %f\n", right.data[0]);
		printf("dim: %d\n", dim.data[0]);
		printf("cellid: %d\n", cellid[100].data[0]);
		printf("position: %f\n", position[100].data[0]);


		hipLaunchKernelGGL(HIP_KERNEL_NAME(ComputeForce<FPTYPE>), dim3(block), dim3(THREADS_PER_BLOCK), 0, 0,
			nAtoms, dt, fmt2v, mass, left, right, dim, cellid, /*locator, */position, v, force);


		hipErrorCheck(hipGetLastError());
		hipErrorCheck(hipDeviceSynchronize());
	}
};

template struct direct_truncation_op<float, device::DEVICE_GPU>;
template struct direct_truncation_op<double, device::DEVICE_GPU>;
}

