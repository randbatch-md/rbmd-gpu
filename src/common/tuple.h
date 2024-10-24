#pragma once

#include <thrust/transform.h>
#include <thrust/device_vector.h>
#include <thrust/tuple.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/functional.h>
#include <utility>  // For std::index_sequence


// 辅助函数：对元组中的元素进行加和
template<typename Tuple, std::size_t... I>
__device__ rbmd::Real sum_tuple(const Tuple& forces_tuple, std::index_sequence<I...>) {
  return (thrust::get<I>(forces_tuple) + ...);  // 逐元素加和
}

// 可变参数模板：对任意数量的力进行加和
template<typename Result, typename... Forces>
void sum_forces_direction(Result& result, Forces&... forces)
{
  // 使用 thrust::transform 并行计算，将多个力相加
  auto zip_begin = thrust::make_zip_iterator(thrust::make_tuple(forces.begin()...));
  auto zip_end = thrust::make_zip_iterator(thrust::make_tuple(forces.end()...));

  thrust::transform(
      zip_begin, zip_end, result.begin(),
      [] __device__ (auto forces_tuple) {
          // 使用索引展开元组并进行加和
          constexpr std::size_t num_forces = thrust::tuple_size<decltype(forces_tuple)>::value;
          return sum_tuple(forces_tuple, std::make_index_sequence<num_forces>{});
      }
  );
}

template<typename... Forces>
void sum_all_forces(
    thrust::device_vector<rbmd::Real>& result_f,
    Forces&... forces)
{
  sum_forces_direction(result_f, forces...);
}
