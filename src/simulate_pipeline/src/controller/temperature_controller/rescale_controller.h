#pragma once
#include "../executioner/executioner.h"
#include "temperature_controller.h"

class RescaleController : public TemperatureController {
 public:
  RescaleController();
  virtual ~RescaleController();

  void Init() override;
  void Update() override;

  /**
   * @brief Calculate the current stage temperature
   */
  void ComputeTemp()override;

  /**
   * @brief Update current speed through temperature
   */
  void UpdataVelocity();

 private:
  rbmd::Id _num_atoms;
  rbmd::Real _mvv2e;
  rbmd::Real _kB;
  rbmd::Real _temp_sum;
  rbmd::Real _temp;
  rbmd::Real* _d_temp_contrib;
};