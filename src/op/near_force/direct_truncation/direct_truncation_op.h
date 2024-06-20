#pragma once
#include "base/device_types.h"

namespace op
{

template <typename FPTYPE, typename Device>
struct direct_truncation_op
{
	void operator()();
};

}