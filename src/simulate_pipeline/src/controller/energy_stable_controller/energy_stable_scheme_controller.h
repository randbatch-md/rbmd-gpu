#pragma once
#include "types.h"
#include "model/device_data.h"
#include "model/md_data.h"
#include "model/structure_info_data.h"


class EnergyStableSchemeController
{
public:
  EnergyStableSchemeController();
  virtual ~EnergyStableSchemeController();

  void Init() ;
  void Update();

  /**
   * @brief Calculate the current stage temperature
   */
  void ComputeTemp();

  /**
   * @brief Update current speed through temperature
   */
  void UpdataVelocity();

private:
  rbmd::Real _dt;
  rbmd::Real _mvv2e;
  rbmd::Real _kB;

  //read
  rbmd::Real _temperature_start;
  rbmd::Real _temperature_stop;
  rbmd::Real _temperature_damp;

  //compute
  rbmd::Real* _d_temp_contrib;
  rbmd::Real  _temp_sum;
  rbmd::Real _temperature;

  std::shared_ptr<StructureInfoData> _structure_info_data;
  std::shared_ptr<DeviceData> _device_data;
};