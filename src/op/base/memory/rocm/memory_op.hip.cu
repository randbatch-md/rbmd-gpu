#include "base/memory/memory_op.h"
#include "base/rocm.h"

namespace op {

	template <typename FPTYPE>
	struct delete_memory_op<FPTYPE, device::DEVICE_GPU>
	{
		void operator()(FPTYPE* arr)
		{
			hipErrorCheck(hipFree(arr));
		}
	};

	template <typename FPTYPE>
	struct resize_memory_op<FPTYPE, device::DEVICE_GPU>
	{
		void operator()(FPTYPE*& arr, const size_t size)
		{
			if (nullptr != arr) {
				delete_memory_op<FPTYPE, device::DEVICE_GPU>(arr);
			}

			hipErrorCheck(hipMalloc((void**)&arr, sizeof(FPTYPE) * size));
		}
	};

	template <typename FPTYPE>
	struct sync_memory_h2d_op<FPTYPE, device::DEVICE_GPU>
	{
		void operator()(FPTYPE* arr_d, const FPTYPE* arr_s, const size_t size)
		{
			hipErrorCheck(hipMemcpy(arr_d, arr_s, sizeof(FPTYPE) * size, hipMemcpyHostToDevice));
		}
	};

	template <typename FPTYPE>
	struct sync_memory_d2h_op<FPTYPE, device::DEVICE_GPU>
	{
		void operator()(FPTYPE* arr_d, const FPTYPE* arr_s, const size_t size)
		{
			hipErrorCheck(hipMemcpy(arr_d, arr_s, sizeof(FPTYPE) * size, hipMemcpyDeviceToHost));
		}
	};

	template struct delete_memory_op <float, device::DEVICE_GPU>;
	template struct delete_memory_op <double, device::DEVICE_GPU>;

	template struct resize_memory_op <float, device::DEVICE_GPU>;
	template struct resize_memory_op <double, device::DEVICE_GPU>;

	template struct sync_memory_h2d_op <float, device::DEVICE_GPU>;
	template struct sync_memory_h2d_op <double, device::DEVICE_GPU>;

	template struct sync_memory_d2h_op <float, device::DEVICE_GPU>;
	template struct sync_memory_d2h_op <double, device::DEVICE_GPU>;
}

