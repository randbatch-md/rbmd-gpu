#include "near_force/direct_truncation/direct_truncation_op.h"
#include <hip/hip_runtime.h>

namespace op
{

#define THREADS_PER_BLOCK 256

__global__ void test_device()
{
	printf("device::test_device()");
}

void test()
{
	int ng = 1;
	int block = (ng + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
	hipLaunchKernelGGL(test_device, dim3(block), dim3(THREADS_PER_BLOCK),0,0);
}

//template<typename FPTYPE>
//void direct_truncation_op<FPTYPE, device::DEVICE_GPU>::operator()()
//{
//	int ng = 1;
//	int block = (ng + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
//	hipLaunchKernelGGL(test<FPTYPE>, dim3(block), dim3(THREADS_PER_BLOCK),0,0);
//}
//
//template struct direct_truncation_op<float, device::DEVICE_GPU>;
//template struct direct_truncation_op<double, device::DEVICE_GPU>;
}

