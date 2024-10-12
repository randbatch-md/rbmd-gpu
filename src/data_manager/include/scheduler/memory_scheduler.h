#pragma once
#include <hip/hip_runtime.h>

#include "common/object.h"

// #include "../data_manager/include/data_manager.h"
#include "../data_manager/include/model/device_data.h"
#include "../data_manager/include/model/md_data.h"

class MemoryScheduler : public Object {
 public:
  /**
   * @brief constructor
   */
  MemoryScheduler();
  virtual ~MemoryScheduler() = default;
  /**
   * @brief async memeory host to device
   * @return error code
   */
  virtual bool asyncMemoryH2D();

  /**
   * @brief async memory device to host
   * @return error code
   */
  virtual bool asyncMemoryD2H();

 protected:
  // DataManager& _data_manager;
  std::shared_ptr<DeviceData>& _device_data;
  const std::shared_ptr<MDData>& _md_data;
  const std::shared_ptr<StructureData>& _structure_data;
  const std::shared_ptr<StructureInfoData>& _structure_info_data;
  const std::shared_ptr<ForceFieldData>& _force_field_data;
};