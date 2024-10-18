#pragma once
#include <memory>

#include "data_manager.h"
#include "model/device_data.h"
#include "model/structure_info_data.h"

class Force {
 public:
  Force() {
    this->_structure_info_data =
        DataManager::getInstance().getMDData()->_structure_info_data;
    this->_device_data = DataManager::getInstance().getDeviceData();
  };
  virtual ~Force() = default;

  // virtual void Update()=0;
  virtual void Init() {};
  virtual void Execute() = 0;

 protected:
  std::shared_ptr<StructureInfoData> _structure_info_data;
  std::shared_ptr<DeviceData> _device_data;
};