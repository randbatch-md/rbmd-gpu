#include "base/device_types.h""
#include "near_force/direct_truncation/direct_truncation_op.h"
#include <iostream>

namespace op
{

template <typename FPTYPE>
__device__ __inline__
void test()
{
	std::cout << "device::test()" << std::endl;
}

template<typename FPTYPE>
void direct_truncation_op<FPTYPE, DEVICE_GPU>::operator()()
{
	test();
}

template struct direct_truncation_op<float, DEVICE_GPU>;
template struct direct_truncation_op<double, DEVICE_GPU>;
}

