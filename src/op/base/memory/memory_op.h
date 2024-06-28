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

	template <typename FPTYPE, typename DEVICE>
	struct sync_memory_h2d_op
	{
		void operator()(FPTYPE* arr_d, const FPTYPE* arr_s, const size_t size);
	};

	template <typename FPTYPE, typename DEVICE>
	struct sync_memory_d2h_op
	{
		void operator()(FPTYPE* arr_d, const FPTYPE* arr_s, const size_t size);
	};
 }