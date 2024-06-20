#pragma once
#include <hip/hip_runtime.h>
#include "base/device_types.h"

namespace op
{

template <typename FPTYPE, typename Device>
struct direct_truncation_op
{
	void operator()();
};

extern template struct direct_truncation_op<float, DEVICE_GPU>;
extern template struct direct_truncation_op<double, DEVICE_GPU>;
}