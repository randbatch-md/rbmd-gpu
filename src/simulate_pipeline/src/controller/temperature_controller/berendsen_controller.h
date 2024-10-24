#pragma once
#include "temperature_controller.h"

class BerendsenController : public TemperatureController {
 public:
  BerendsenController();
  virtual ~BerendsenController();

  void Init() override;
  void Update() override;

  /**
   * @brief Calculate the current stage temperature
   */
  void ComputeTemp() override;

  /**
   * @brief Update current speed through temperature
   */
  void UpdataVelocity();

 private:
  rbmd::Id _num_atoms;
  rbmd::Real _dt;
  rbmd::Real _mvv2e;
  rbmd::Real _kB;
  rbmd::Real _temp_sum;
  rbmd::Real _temp;

  rbmd::Real _Tdamp;
  rbmd::Real* _d_temp_contrib;
};