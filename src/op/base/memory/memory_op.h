#pragma once
#include "base/device_types.h"
#include "types.h"

namespace op {

	template <typename FPTYPE, typename DEVICE>
	struct delete_memory_op
	{
		void operator()(FPTYPE* arr);
	};

	template <typename FPTYPE, typename DEVICE>
	struct resize_memory_op
	{
		void operator()(FPTYPE*& arr, const size_t size);
	};

	template <typename FPTYPE, typename DEVICE_S, typename DEVICE_D>
	struct sync_memory_op
	{
		void operator()(const FPTYPE* arr_s, FPTYPE arr_d, const size_t size);
	};
 }