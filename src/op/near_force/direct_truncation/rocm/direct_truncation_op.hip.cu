#include "near_force/direct_truncation/direct_truncation_op.h"
#include <hip/hip_runtime.h>

namespace op
{

#define THREADS_PER_BLOCK 256

template<typename FPTYPE>
__global__ void test_LJ()
{
	printf("device::test_device()");
}

void LJ()
{
	int ng = 1;
	int block = (ng + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
	hipLaunchKernelGGL(test_LJ<FPTYPE>, dim3(block), dim3(THREADS_PER_BLOCK),0,0);
}

//template<typename FPTYPE>
//void direct_truncation_op<FPTYPE, device::DEVICE_GPU>::operator()()
//{
//	int ng = 1;
//	int block = (ng + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
//	hipLaunchKernelGGL(test_LJ<FPTYPE>, dim3(block), dim3(THREADS_PER_BLOCK),0,0);
//}
//
//template struct direct_truncation_op<float, device::DEVICE_GPU>;
//template struct direct_truncation_op<double, device::DEVICE_GPU>;
}

