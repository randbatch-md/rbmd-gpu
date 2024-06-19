#pragma once
#include "base/device_types.h"

namespace op
{

template <typename FPTYPE, DeviceType Device>
struct direct_truncation_op
{
	void operator()();
};

extern template struct direct_truncation_op<float, DEVICE_GPU>;
extern template struct direct_truncation_op<double, DEVICE_GPU>;
}