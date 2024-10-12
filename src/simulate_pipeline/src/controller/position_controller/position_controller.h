#pragma once
#include "data_manager.h"
#include "model/device_data.h"
#include "model/md_data.h"
#include "model/structure_info_data.h"
#include "shake_controller.h"
class PositionController {
 public:
  PositionController()
      : _device_data(DataManager::getInstance().getDeviceData()),
        _structure_info_data(
            DataManager::getInstance().getMDData()->_structure_info_data){};

  virtual ~PositionController() = default;

  /**
   * @brief Update Position
   */
  virtual void Update() = 0;

  /**
   * @brief Parameters and objects required for initializing the Position
   * controller
   */
  virtual void Init() = 0;

 protected:
  std::shared_ptr<StructureInfoData> _structure_info_data;
  std::shared_ptr<ShakeController> _shake_controller;
  std::shared_ptr<DeviceData> _device_data;
};