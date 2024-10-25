#pragma once
#include <memory>

#include "data_manager.h"
#include "device_types.h"
#include "model/device_data.h"
#include "model/md_data.h"
#include "model/structure_info_data.h"
#include "types.h"
class TemperatureController {
 public:
  TemperatureController()
      : _device_data(DataManager::getInstance().getDeviceData()),
        _structure_info_data(
            DataManager::getInstance().getMDData()->_structure_info_data){};

  virtual ~TemperatureController() = default;

  /**
   * @brief Update Temperature
   */
  virtual void Update() = 0;

  virtual void ComputeTemp() {};
  /**
   * @brief Parameters and objects required for initializing the temperature
   * controller
   */
  virtual void Init() = 0;

 protected:
  std::shared_ptr<StructureInfoData> _structure_info_data;
  std::shared_ptr<DeviceData> _device_data;
};