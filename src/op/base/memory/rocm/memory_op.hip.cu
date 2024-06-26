#include "base/memory/memory_op.h"
#include "base/rocm.h"
#include "base/device_types.h"
#include "model/types.h"

namespace op {

	template <typename FPTYPE, device::DEVICE_GPU>
	struct delete_memory_op
	{
		void operator()(FPTYPE* arr) 
		{
			hipErrorCheck(hipFree(arr));
		}
	};

	template <typename FPTYPE, device::DEVICE_GPU>
	struct resize_memory_op
	{
		void operator()(FPTYPE*& arr, const size_t size)
		{
			if (nullptr != arr) {
				delete_memory_op<FPTYPE, device::DEVICE_GPU>(arr);
			}

			hipErrorCheck(hipMalloc((void**)&arr, sizeof(FPTYPE) * size));
		}
	};

	template <typename FPTYPE, device::DEVICE_CPU, device::DEVICE_GPU>
	struct sync_memory_op
	{
		void operator()(const FPTYPE* arr_s, FPTYPE arr_d, const size_t size)
		{
			hipErrorCheck(hipMemcpy(arr_s, arr_d, sizeof(FPTYPE) * size, hipMemcpyHostToDevice));
		}
	};

	template <typename FPTYPE, device::DEVICE_GPU, device::DEVICE_CPU>
	struct sync_memory_op
	{
		void operator()(const FPTYPE* arr_s, FPTYPE arr_d, const size_t size)
		{
			hipErrorCheck(hipMemcpy(arr_s, arr_d, sizeof(FPTYPE) * size, hipMemcpyDeviceToHose));
		}
	};

	template struct delete_memory_op <float, device::DEVICE_GPU>;
	template struct delete_memory_op <double, device::DEVICE_GPU>;

	template struct resize_memory_op <float, device::DEVICE_GPU>;
	template struct resize_memory_op <double, device::DEVICE_GPU>;

	template struct sync_memory_op <float, device::DEVICE_CPU, device::DEVICE_GPU>;
	template struct sync_memory_op <double, device::DEVICE_CPU, device::DEVICE_GPU>;

	template struct sync_memory_op <float, device::DEVICE_GPU, device::DEVICE_CPU>;
	template struct sync_memory_op <double, device::DEVICE_GPU, device::DEVICE_CPU>;
}

