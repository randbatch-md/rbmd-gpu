#pragma once
#include "../common/object.h"
#include "force.h"
#include "position_controller.h"
#include "shake_controller.h"
#include "temperature_controller.h"
#include "pressure_controller.h"
#include "nose_hoover_controller.h"
#include "velocity_controller.h"
#include "energy_stable_scheme_controller.h"

class Ensemble : public Object {
 public:
  Ensemble(){};
  virtual ~Ensemble() = default;

  /**
   *@brief init current ensemble
   */
  virtual void Init() = 0;

  /**
   *@brief Preprocessing settings before solving
   */
  virtual void Presolve() = 0;

  /**
   * @brief Specific execution calculation part
   */
  virtual void Solve() = 0;

  /**
   * @brief Post processing after the completion of the current calculation step
   */
  virtual void Postsolve() = 0;

  /**
   * @brief The overall calculation process of the current calculation step
   * @return Determine if there is an error here. If executed normally, return
   * 0. If not, return other values
   */
  int Run() {
    Presolve();
    Solve();
    Postsolve();

    return 0;
  };

 protected:
  std::shared_ptr<PositionController> _position_controller;
  std::shared_ptr<VelocityController> _velocity_controller;
  std::shared_ptr<Force> _force_controller;
  std::shared_ptr<ShakeController> _shake_controller;
  std::shared_ptr<TemperatureController> _temperature_controller;
  std::shared_ptr<PressureController> _pressure_controller;
  std::shared_ptr<NoseHooverController> _NoseHoover_controller;
  std::shared_ptr<EnergyStableSchemeController> _energy_stable_scheme_controller;

  std::string _temp_ctrl_type, _press_ctrl_type;
  std::string  _integration_type;
  size_t _total_steps = 0;            // 累计总步数
  double _total_time_vv = 0.0;        // VV算法总耗时
  double _total_time_beeman = 0.0;    // Beeman算法总耗时
  double _total_time_prk3c = 0.0;     // PRK3C算法总耗时
};