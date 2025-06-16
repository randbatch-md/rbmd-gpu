#pragma once
#include "../../common/device_types.h"
#include "../../common/types.h"
#include "../../data_manager/include/model/box.h"

namespace op {

//

template <typename DEVICE>
struct FixRBLForceSingleOp
{
  void operator()(
          const rbmd::Id num_atoms,
          const rbmd::Real corr_value_single,
          rbmd::Real* f_single);
};

// template <typename DEVICE>
// struct FixRBLForceOp
// {
//   void operator()(
//           const rbmd::Id num_atoms,
//           const rbmd::Real corr_value_x,
//           const rbmd::Real corr_value_y,
//           const rbmd::Real corr_value_z,
//           rbmd::Real* fx,
//           rbmd::Real* fy,
//           rbmd::Real* fz);
// };

template <>
struct FixRBLForceSingleOp<device::DEVICE_GPU>
{
  void operator()(
          const rbmd::Id num_atoms,
          const rbmd::Real corr_value_single,
          rbmd::Real* f_single);
};

// template <>
// struct FixRBLForceOp<device::DEVICE_GPU>
// {
//   void operator()(
//           const rbmd::Id num_atoms,
//           const rbmd::Real corr_value_x,
//           const rbmd::Real corr_value_y,
//           const rbmd::Real corr_value_z,
//           rbmd::Real* fx,
//           rbmd::Real* fy,
//           rbmd::Real* fz);
// };

}