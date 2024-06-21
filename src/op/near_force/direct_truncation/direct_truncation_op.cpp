#include "direct_truncation_op.h"
#include "device_types.h"
#include<iostream>

namespace op
{
	template<typename FPTYPE>
	void LJ()
	{
		std::cout << "CPU LJ" << std::endl;
	}
	template void LJ<float>();
	template void LJ<double>();

	//template <typename FPTYPE>
	//struct direct_truncation_op<FPTYPE, device::DEVICE_CPU>
	//{
	//	void operator()()
	//	{
	//		std::cout << "CPU" << std::endl;
	//	}
	//};


	//template struct direct_truncation_op<float, device::DEVICE_CPU>;
	//template struct direct_truncation_op<double, device::DEVICE_CPU>;
}