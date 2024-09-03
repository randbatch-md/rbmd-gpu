#pragma once
#include "box.h"
#include "common/types.h"

namespace op {
template <typename DEVICE>
struct InitEndIndexOp {
  void operator()(rbmd::Id* neighbor_num, rbmd::Id* start_index,
                  rbmd::Id* end_index, rbmd::Id total_atom_num);
};
} // namespace op
