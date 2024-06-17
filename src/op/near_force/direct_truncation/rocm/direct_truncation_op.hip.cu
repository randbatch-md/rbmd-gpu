#include "base/device_types.h""
#include "near_force/direct_truncation/direct_truncation_op.h"
#include <iostream>

namespace op
{

template<typename FPTYPE>
struct direct_truncation_op<FPTYPE, DEVICE_GPU>
{
	void operator()()
	{
		std::cout << "GPU" << std::endl;
	}
};

template struct direct_truncation_op<float, DEVICE_GPU>;
template struct direct_truncation_op<double, DEVICE_GPU>;
}

