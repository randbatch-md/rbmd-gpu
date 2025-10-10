#pragma once
#include <memory>

#include "data_manager.h"
#include "device_types.h"
#include "model/device_data.h"
#include "model/md_data.h"
#include "model/structure_info_data.h"
#include "types.h"
extern int test_current_step;

class TemperatureController {
 public:
  TemperatureController()
      : _device_data(DataManager::getInstance().getDeviceData())
      ,_structure_info_data(DataManager::getInstance().getMDData()->_structure_info_data)
      , _temp_sum(0){};

  virtual ~TemperatureController() = default;

  /**
   * @brief Update Temperature
   */
  virtual void Update() = 0;

  virtual void ComputeTemperature() {};
  /**
   * @brief Parameters and objects required for initializing the temperature
   * controller
   */
  virtual void Init() = 0;

  void ComputeTempTargetInit() {
    if (_temperature_stop == _temperature_start) //Thermostatic simulation
    {
      _t_target = _temperature_stop= _temperature_start;
    }
    else                    //anisothermal simulation
    {
      auto currentstep = test_current_step;
      auto beginstep = 0;
      auto endstep = DataManager::getInstance().getConfigData()->
        Get<rbmd::Real>("num_steps", "execution");

      rbmd::Real delta = currentstep - beginstep;

      if (delta != 0.0)
      {
        delta = delta / static_cast<rbmd::Real>(endstep - beginstep);
      }

      _t_target = _temperature_start + delta * (_temperature_stop - _temperature_start);
    }
  }

 protected:
  std::shared_ptr<StructureInfoData> _structure_info_data;
  std::shared_ptr<DeviceData> _device_data;

  rbmd::Real _temp_sum;
  rbmd::Real _temperature;

  rbmd::Real _temperature_start =0.0;   //read
  rbmd::Real _temperature_stop =0.0;    //read
  rbmd::Real _temperature_damp =0.0;    //read

  rbmd::Real  _t_target;        //compute
  bool _com_bias = false;

  std::string _group_name = "all";
  Real3 _vbias;


};