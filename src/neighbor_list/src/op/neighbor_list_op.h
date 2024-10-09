#pragma once
#include "common/device_types.h"
#include "common/types.h"
namespace op {
template <typename DEVICE>
struct InitEndIndexOp {
  void operator()(rbmd::Id* neighbor_num, rbmd::Id* start_index,
                  rbmd::Id* end_index, rbmd::Id total_atom_num);
};

template <>
struct InitEndIndexOp<device::DEVICE_GPU> {
  void operator()(rbmd::Id* neighbor_num, rbmd::Id* start_index,
                  rbmd::Id* end_index, rbmd::Id total_atom_num);
};
}  // namespace op
