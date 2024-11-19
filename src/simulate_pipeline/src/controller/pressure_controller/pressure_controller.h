#pragma once
#include <memory>

#include "data_manager.h"
#include "device_types.h"
#include "model/device_data.h"
#include "model/md_data.h"
#include "model/structure_info_data.h"
#include "types.h"
class PressureController {
public:
  PressureController()
      : _structure_info_data(
           DataManager::getInstance().getMDData()->_structure_info_data)
      , _device_data(DataManager::getInstance().getDeviceData())
      ,
       _box(DataManager::getInstance().getMDData()->_box)
      , _pressure(0){};

  virtual ~PressureController() = default;

  /**
   * @brief Update Temperature
   */
  virtual void Update() = 0;

  virtual void ComputePressure() {};
  /**
   * @brief Parameters and objects required for initializing the pressure
   * controller
   */
  virtual void Init() = 0;

protected:
  std::shared_ptr<StructureInfoData> _structure_info_data;
  std::shared_ptr<DeviceData> _device_data;
  std::shared_ptr<Box>  _box;

  rbmd::Real _pressure;
};