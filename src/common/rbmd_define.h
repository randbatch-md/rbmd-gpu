#pragma once
#if defined(__CUDA)
    #include <cuda_runtime.h>
    #include <cuda_runtime_api.h>
    #include <thrust/device_vector.h>
    #include <cub/cub.cuh>
    #include<cmath>
    #include <thrust/gather.h>
    #include <thrust/sort.h>
#elif defined (__ROCM)
    #include <hip/hip_runtime.h>
    #include <hip/hip_runtime_api.h>
    #include <thrust/device_vector.h>
    #include <hipcub/hipcub.hpp>
    #include <thrust/gather.h>
    #include <thrust/sort.h>
    #include <hipcub/backend/rocprim/device/device_radix_sort.hpp>
    #include <hipcub/backend/rocprim/iterator/counting_input_iterator.hpp>
#else
     #error "This code must be compiled with either HIP or CUDA."
#endif


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
#define RBMD_TRUE (1)
#define RBMD_FALSE (0)

#ifndef TARGET_DCU
#define MIN_NBNUM \
(128)  /// CUDA AMD6800xt 96 DCU 128   TODO kernel us it  can use warpSize?
#define WARP_SIZE (32)  /// CUDA AMD6800xt 32  DCU 64   TODO
#else
#define MIN_NBNUM (128)
#define WARP_SIZE (64)
#endif

#if USE_DOUBLE
    #if defined(__ROCM)
        #define REAL_DATA(vec) (vec.data)
    #elif defined(__CUDA)
        #define REAL_DATA(vec) (const_cast<double*>(reinterpret_cast<const double*>(&vec)))
    #endif
typedef double3 Real3;
typedef double2 Real2;
#define make_Real3 make_double3
#define make_Real2 make_double2
#define POW pow
#define CEIL ceil
#define FLOOR floor
#define SQRT sqrt
#define ERF erf
#define EXP exp
#define COS cos
#define SIN sin
#define ROUND round
#define LOG log
#define ABS fabs
#define ACOS acos
#else
    #if defined(__ROCM)
        #define REAL_DATA(vec) (vec.data)
    #elif defined(__CUDA)
        #define REAL_DATA(vec) (const_cast<float*>(reinterpret_cast<const float*>(&vec)))
    #endif
typedef float3 Real3;
typedef float2 Real2;
#define make_Real3 make_float3
#define make_Real2 make_float2
#define POW powf
#define CEIL ceilf
#define FLOOR floorf
#define SQRT sqrtf
#define ERF erff
#define EXP expf
#define COS cosf
#define SIN sinf
#define ROUND roundf
#define LOG logf
#define ABS fabsf
#define ACOS acosf

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

#if defined(__GNUC__) || defined(__CUDA)  // GCC
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

#if defined(__CUDA)  // NVCC   //TODO   待验证
#define ALIGN(n) __align__(n)
#elif defined(__GNUC__)  // GCC
#define ALIGN(n) __attribute__((aligned(n)))
#elif defined(_MSC_VER)  // MSVC
#define ALIGN(n) __declspec(align(n))
#else
#error "Please provide a definition for ALIGN macro for your host compiler!"
#endif

#if defined (__CUDA)
    #define MALLOC cudaMalloc
    #define MALLOCHOST(ptr, size) cudaHostAlloc((void**)ptr, size,cudaHostAllocDefault)
    #define MEMCPY cudaMemcpy
    #define H2D cudaMemcpyHostToDevice
    #define H2H cudaMemcpyHostToHost
    #define D2H cudaMemcpyDeviceToHost
    #define D2D cudaMemcpyDeviceToDevice
    #define FREE cudaFree
    #define MEMSET cudaMemset
    #define REDUCE cub::DeviceReduce::Sum
    #define ERROR_T cudaError_t
    #define SUCCESS cudaSuccess
    #define GETERRORSTRING cudaGetErrorString
    #define GETERRORNAME cudaGetErrorName
    #define LASTERROR cudaPeekAtLastError
    #define EXCLUSIVESUM cub::DeviceScan::ExclusiveSum
    #define WARPREDUCE cub::WarpReduce
    #define WARPSCAN cub::WarpScan
    #define SHUFFLEINDEX cub::ShuffleIndex
    #define BLOCKREDUCE cub::BlockReduce
#elif defined (__ROCM)
    #define MALLOC hipMalloc
    #define MALLOCHOST hipHostMalloc
    #define MEMCPY hipMemcpy
    #define H2D hipMemcpyHostToDevice
    #define H2H hipMemcpyHostToHost
    #define D2H hipMemcpyDeviceToHost
    #define D2D hipMemcpyDeviceToDevice
    #define FREE hipFree
    #define MEMSET hipMemset
    #define REDUCE hipcub::DeviceReduce::Sum
    #define ERROR_T hipError_t
    #define SUCCESS hipSuccess
    #define GETERRORSTRING hipGetErrorString
    #define GETERRORNAME cudaGetErrorName
    #define LASTERROR hipPeekAtLastError
    #define EXCLUSIVESUM hipcub::DeviceScan::ExclusiveSum
    #define WARPREDUCE hipcub::WarpReduce
    #define WARPSCAN hipcub::WarpScan
    #define SHUFFLEINDEX hipcub::ShuffleIndex
    #define BLOCKREDUCE hipcub::BlockReduce
#endif


template <typename T>
static T *raw_ptr(thrust::device_vector<T> &vec) {
  return thrust::raw_pointer_cast(vec.data());
}


#define CHECK_RUNTIME(call) CheckRunTime(call, #call, __LINE__, __FILE__)
static bool CheckRunTime(ERROR_T e, const char* call, int line,
    const char* file)
{
    if (e != SUCCESS) {
        printf("Runtime error %s # %s, code = %s [ %d ] in file %s:%d", call,
            GETERRORSTRING(e), GETERRORNAME(e), e, file, line);
        return false;
    }
    return true;
}

#define CHECK_KERNEL(...)                                               \
  __VA_ARGS__;                                                          \
  do {                                                                  \
    ERROR_T err = LASTERROR();                       \
    if (err != SUCCESS) {                                     \
      printf("Launch Kernel Failed:  %s:%d '%s'\n", __FILE__, __LINE__, \
             GETERRORSTRING(err));                            \
      exit(EXIT_FAILURE);                                               \
    }                                                                   \
  } while (0);


template <typename T>
// d_src_array input array  d_dst outputnum size：input array size
static void ReductionSum(T* d_src_array, T* d_dst, rbmd::Id size) {
    //void* temp = nullptr;
    //size_t temp_bytes = 0;
    //CHECK_RUNTIME(REDUCE(temp, temp_bytes, d_src_array, d_dst,
    //    static_cast<int>(size)));
    //CHECK_RUNTIME(MALLOC(&temp, temp_bytes));
    //CHECK_RUNTIME(REDUCE(temp, temp_bytes, d_src_array, d_dst,
    //    static_cast<int>(size)));
    //CHECK_RUNTIME(FREE(temp));
}


template<typename Tuple, std::size_t... I>
__device__ rbmd::Real SumTuple(const Tuple& forces_tuple, std::index_sequence<I...>) {
  return (thrust::get<I>(forces_tuple) + ...);
}


template<typename Result, typename... Forces>
void SumforcesDirection(Result& result, Forces&... forces)
{

  auto zip_begin = thrust::make_zip_iterator(thrust::make_tuple(forces.begin()...));
  auto zip_end = thrust::make_zip_iterator(thrust::make_tuple(forces.end()...));

  thrust::transform(
      zip_begin, zip_end, result.begin(),
      [] __device__ (auto forces_tuple) {
          constexpr std::size_t num_forces = thrust::tuple_size<decltype(forces_tuple)>::value;
          return SumTuple(forces_tuple, std::make_index_sequence<num_forces>{});
      }
  );
}

template<typename... Forces>
void TransformForces(
    thrust::device_vector<rbmd::Real>& result_f,
    Forces&... forces)
{
  SumforcesDirection(result_f, forces...);
}