#pragma once
#include "base/device_types.h"

namespace op
{
	template<typename FPTYPE>
	void LJ();

template <typename FPTYPE, typename DEVICE>
struct direct_truncation_op
{
	void operator()(int test);
};

}