#pragma once
#include <memory>
#include "data_manager.h"
#include "model/md_data.h"

class Output {
public:
  Output() {
    this->_device_data = DataManager::getInstance().getDeviceData();
    this->_structure_info_data = DataManager::getInstance().getMDData()->_structure_info_data;
  };

  virtual ~Output() = default;

  /**
   * @brief Execute output
   */
  virtual void Execute(int current_timestep) = 0;
  /**
   * @brief Parameters and objects required for initializing the output
   * controller
   */
  virtual void  Init(size_t num_atoms, rbmd::Id interval) = 0;

protected:
  std::shared_ptr<StructureInfoData> _structure_info_data;
  std::shared_ptr<DeviceData> _device_data;
};