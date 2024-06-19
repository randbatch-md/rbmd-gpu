#include "base/device_types.h""
#include "near_force/direct_truncation/direct_truncation_op.h"
#include <iostream>
#include <hip/hip_runtime.h>

namespace op
{

#define THREADS_PER_BLOCK 256

template <typename FPTYPE>
__global__
void test()
{
	std::cout << "device::test()" << std::endl;
}

template<typename FPTYPE>
void direct_truncation_op<FPTYPE, DEVICE_GPU>::operator()()
{
	int ng = 1;
	int block = (ng + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
	hipLaunchKernelGGL(test<FPTYPE>, dim3(block), dim3(THREADS_PER_BLOCK),0,0);
}

template struct direct_truncation_op<float, DEVICE_GPU>;
template struct direct_truncation_op<double, DEVICE_GPU>;
}

