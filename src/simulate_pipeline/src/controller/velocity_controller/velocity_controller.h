#pragma once
// #include "model/device_data.h"
// #include "model/structure_info_data.h"
#include <memory>

#include "data_manager.h"
#include "model/md_data.h"

extern int test_current_step;

class VelocityController {
 public:
  VelocityController() {
    this->_device_data =
        DataManager::getInstance()
            .getDeviceData();  // todo 最开始的初始化待定
                               // _device_data(std::make_shared<DeviceData>())
    this->_structure_info_data =
        DataManager::getInstance()
            .getMDData()
            ->_structure_info_data;  //  todo 最开始的初始化待定
                                     //  _structure_info_data(std::make_shared<StructureInfoData>())
  };

  virtual ~VelocityController() = default;

  /**
   * @brief Update Speed
   */
  virtual void Update() = 0;
  /**
   * @brief Parameters and objects required for initializing the speed
   * controller
   */
  virtual void Init() = 0;
  virtual void Update1() = 0;
  virtual void Update2() = 0;
  virtual void Update3() = 0;
  virtual void Update4() = 0;
  virtual void Updatevl() = 0;
  virtual void Updatebm() = 0;

 protected:
  std::shared_ptr<StructureInfoData> _structure_info_data;
  std::shared_ptr<DeviceData> _device_data;
};