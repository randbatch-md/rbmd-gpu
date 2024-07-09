#pragma once
#include "base/device_types.h"
#include "types.h"

namespace op
{

template <typename FPTYPE, typename DEVICE>
struct direct_truncation_op
{
	void operator()(
		const int& nSteps,
		const int& nAtoms,
		const FPTYPE* dt,
		const FPTYPE* fmt2v,
		const FPTYPE* mass,
		rbmd::Real3* position,
		rbmd::Real3* v,
		rbmd::Real3* force);
};

}