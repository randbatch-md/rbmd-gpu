#pragma once
#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>
#include <thrust/device_vector.h>

#include <hipcub/hipcub.hpp>

#include "types.h"
#if USE_DOUBLE
#define EPSILON 0.0001
#else
#define EPSILON 0.0001f
#endif

#define SAFE_ZONE (1.2)
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define BLOCK_SIZE (256)
#define MAX_GPU_STREAMS (6)
#ifdef AMD_CUDA
#define MIN_NBNUM \
  (96)  /// CUDA AMD6800xt 96 DCU 128   TODO kernel us it  can use warpSize?
#define WARP_SIZE (32)  /// CUDA AMD6800xt 32  DCU 64   TODO
#else
#define MIN_NBNUM (128)
#define WARP_SIZE (64)
#endif
#if USE_DOUBLE
typedef double3 Real3;
#define make_Real3 make_double3
#define POW pow
#define CEIL ceil
#define FLOOR floor
#define SQRT sqrt
#define ERF erf
#define EXP exp
#else
typedef float3 Real3;
#define make_Real3 make_float3
#define POW powf
#define CEIL ceilf
#define FLOOR floorf
#define SQRT sqrtf
#define ERF erff
#define EXP expf
#endif

#if USE_64BIT_IDS
typedef longlong3 Int3;
#define make_Int3 make_longlong3
#else
typedef int3 Int3;
#define make_Int3 make_int3
#endif

// 返回数组需要对齐的大小，n为数组的长度
#define ALIGN_SIZE(type, n) \
  ((sizeof(type) > 4) ? NEXT_POWER_OF_TWO(n) * 8 : NEXT_POWER_OF_TWO(n) * 4)

#if defined(__GNUC__)  // GCC
#define IS_POWER_OF_TWO(x) (((x) & ((x) - 1)) == 0)
#define NEXT_POWER_OF_TWO(n)      \
  ((n) == 0 ? 1                   \
            : (IS_POWER_OF_TWO(n) \
                   ? (n)          \
                   : (1 << (sizeof(n) * 8 - __builtin_clz((n) - 1)))))

#elif defined(_MSC_VER)  // MSVC   TODO： 待验证
#include <intrin.h>
#define IS_POWER_OF_TWO(x) (((x) & ((x) - 1)) == 0)
#define NEXT_POWER_OF_TWO(n)       \
  ((n) == 0                        \
       ? 1                         \
       : (IS_POWER_OF_TWO(n) ? (n) \
                             : (1 << (sizeof(n) * 8 - _lzcnt_u32((n) - 1)))))

#else
#error "Unsupported compiler"
#endif

#if defined(__CUDACC__)  // NVCC   //TODO   待验证
#define ALIGN(n) __align__(n)
#elif defined(__GNUC__)  // GCC
#define ALIGN(n) __attribute__((aligned(n)))
#elif defined(_MSC_VER)  // MSVC
#define ALIGN(n) __declspec(align(n))
#else
#error "Please provide a definition for ALIGN macro for your host compiler!"
#endif

#define MALLOC hipMalloc
#define MALLOCHOST hipHostMalloc
#define MEMCPY hipMemcpy
#define H2D hipMemcpyHostToDevice
#define H2H hipMemcpyHostToHost
#define D2H hipMemcpyDeviceToHost
#define D2D hipMemcpyDeviceToDevice
#define FREE hipFree
#define MEMSET hipMemset

template <typename T>
static T *raw_ptr(thrust::device_vector<T> &vec) {
  return thrust::raw_pointer_cast(vec.data());
}

#define CHECK_RUNTIME(call) CheckHipRuntime(call, #call, __LINE__, __FILE__)

static bool CheckHipRuntime(hipError_t e, const char *call, int line,
                            const char *file) {
  if (e != hipSuccess) {
    printf("CUDA Runtime error %s # %s, code = %s [ %d ] in file %s:%d", call,
           hipGetErrorString(e), hipGetErrorName(e), e, file, line);
    return false;
  }
  return true;
}

#define CHECK_KERNEL(...)                                               \
  __VA_ARGS__;                                                          \
  do {                                                                  \
    hipError_t hip_status = hipPeekAtLastError();                       \
    if (hip_status != hipSuccess) {                                     \
      printf("Launch Kernel Failed:  %s:%d '%s'\n", __FILE__, __LINE__, \
             hipGetErrorString(hip_status));                            \
      exit(EXIT_FAILURE);                                               \
    }                                                                   \
  } while (0);

template <typename T>
// d_src_array input array  d_dst outputnum size：input array size
static void ReductionSum(T *d_src_array, T *d_dst, rbmd::Id size) {
  void *temp = nullptr;
  size_t temp_bytes = 0;
  CHECK_RUNTIME(hipcub::DeviceReduce::Sum(temp, temp_bytes, d_src_array, d_dst,
                                          static_cast<int>(size)));
  CHECK_RUNTIME(MALLOC(&temp, temp_bytes));
  CHECK_RUNTIME(hipcub::DeviceReduce::Sum(temp, temp_bytes, d_src_array, d_dst,
                                          static_cast<int>(size)));
  CHECK_RUNTIME(FREE(temp));
}
