#pragma once
#include "types.h"
#include "vector"
#include "data_manager.h"
#include "model/md_data.h"

class ShakeController {
 public:
  ShakeController();
  virtual ~ShakeController() = default;

  void Init();
  void ShakeA();
  void ShakeB();

protected:
  std::shared_ptr<StructureInfoData> _structure_info_data;
  std::shared_ptr<DeviceData> _device_data;

  rbmd::Id _num_angle;
  rbmd::Real _dt;
  rbmd::Real _fmt2v;

  std::vector<rbmd::Real> _shake_px;
  std::vector<rbmd::Real> _shake_py;
  std::vector<rbmd::Real> _shake_pz;

  std::vector<rbmd::Real> _shake_vx;
  std::vector<rbmd::Real> _shake_vy;
  std::vector<rbmd::Real> _shake_vz;
};