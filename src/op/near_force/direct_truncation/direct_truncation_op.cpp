#include "direct_truncation_op.h"
#include "device_types.h"
#include<iostream>

namespace op
{
	template <typename FPTYPE>
	struct direct_truncation_op<FPTYPE, device::DEVICE_CPU>
	{
		//to do openmp
		void operator()(
			const int& nSteps,
			const int& nAtoms,
			const FPTYPE* dt,
			const FPTYPE* fmt2v,
			const FPTYPE* mass,
			rbmd::Real3* d_position,
			rbmd::Real3* v,
			rbmd::Real3* force)
		{
			std::cout << "cpu: to do openmp!" << std::endl;
		}
	};


	template struct direct_truncation_op<float, device::DEVICE_CPU>;
	template struct direct_truncation_op<double, device::DEVICE_CPU>;
}