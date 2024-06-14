#include "direct_truncation_op.h"
#include "base/device_types.h"

#include<iostream>

namespace op
{

	template <typename FPTYPE>
	struct direct_truncation_op<FPTYPE, DEVICE_CPU>
	{
		void operator()()
		{
			std::cout << "CPU" << std::endl;
		}
	};


	template struct direct_truncation_op<float, DEVICE_CPU>;
	template struct direct_truncation_op<double, DEVICE_CPU>;
}