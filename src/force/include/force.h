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
    this->_box = DataManager::getInstance().getMDData()->_box;
  };
  virtual ~Force() = default;

  // virtual void Update()=0;
  virtual void Init() {};
  virtual void Execute() = 0;
  virtual void EvaluatePotentialenergy(){};


 protected:
  std::shared_ptr<StructureInfoData> _structure_info_data;
  std::shared_ptr<DeviceData> _device_data;

  std::shared_ptr<Box> _box;
};