#pragma once
#include "device_types.h"
#include "types.h"
#include "model/box.h"
namespace op {

    template <typename DEVICE>
    struct ShakeAOp {
      void operator()(const rbmd::Id num_angle, 
                      const rbmd::Real dt,
                      const rbmd::Real fmt2v, 
                      Box* box,
                      const rbmd::Real min_x, const rbmd::Real min_y,
                      const rbmd::Real min_z, const rbmd::Real max_x,
                      const rbmd::Real max_y, const rbmd::Real max_z,
                      const rbmd::Real* mass, 
                      const rbmd::Id* atoms_type,
                      const Id3* angle_id_vec,
                      rbmd::Real* shake_px,
                      rbmd::Real* shake_py,
                      rbmd::Real* shake_pz,
                      rbmd::Real* shake_vx,
                      rbmd::Real* shake_vy,
                      rbmd::Real* shake_vz,
                      const rbmd::Real* fx,
                      const rbmd::Real* fy, 
                      const rbmd::Real* fz, 
                      rbmd::Real* flag_px,
                      rbmd::Real* flag_py, 
                      rbmd::Real* flag_pz);
    };
    
    template <>
    struct ShakeAOp<device::DEVICE_GPU> {
      void operator()(const rbmd::Id num_angle,
                      const rbmd::Real dt,
                      const rbmd::Real fmt2v,
                      Box* box,
                      const rbmd::Real min_x, const rbmd::Real min_y,
                      const rbmd::Real min_z, const rbmd::Real max_x,
                      const rbmd::Real max_y, const rbmd::Real max_z,
                      const rbmd::Real* mass,
                      const rbmd::Id* atoms_type,
                      const Id3* angle_id_vec,
                      rbmd::Real* shake_px,
                      rbmd::Real* shake_py,
                      rbmd::Real* shake_pz,
                      rbmd::Real* shake_vx,
                      rbmd::Real* shake_vy,
                      rbmd::Real* shake_vz,
                      const rbmd::Real* fx,
                      const rbmd::Real* fy,
                      const rbmd::Real* fz,
                      rbmd::Real* flag_px,
                      rbmd::Real* flag_py,
                      rbmd::Real* flag_pz);
    };
}  // namespace op