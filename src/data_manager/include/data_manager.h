#pragma once
#include <memory>

#include "config_data.h"
//#include "./scheduler/memory_scheduler.h"
// #include "../data_manager/include/model/post_process_data.h"
#include "../data_manager/include/model/device_data.h"
#include "../data_manager/include/model/md_data.h"
#include "../data_manager/include/scheduler/memory_scheduler.h"

class MDData;
class PostProcessData;
// class ConfigData;
// class DeviceData;
//class MemoryScheduler;

class DataManager {
 public:
  /**
   * @brief The copy constructor and assignment operator are not open to the
   * public
   * @param
   */
  DataManager(const DataManager&) = delete;

  DataManager& operator=(const DataManager&) = delete;

 private:
  /**
   * @brief Constructors and destructors are not open to the public
   */
  DataManager()
      : _config_data(std::make_shared<ConfigData>("rbmd.json")),
        _md_data(std::make_shared<MDData>()),
        _device_data(std::make_shared<DeviceData>())
        //, _memory_scheduler(std::make_shared<MemoryScheduler>())
        {};  // �����ʼ����_config_data ԭ����DataManager() = default;

  ~DataManager() = default;

 public:
  /**
   * @brief get Singleton
   * @return static data_manager
   */
  static DataManager& getInstance() {
    static DataManager data_manager;
    return data_manager;
  }

  /**
   * @brief get _md_data
   * @return _md_data
   */
  auto& getMDData() const { return _md_data; }

  /**
   * @brief get _postProcess_data
   * @return _postProcess_data
   */
  auto& getPostProcessData() const { return _postProcess_data; }

  /**
   * @brief get _config_data
   * @return _config_data
   */
  auto& getConfigData() const { return _config_data; }

  /**
   * @brief get _device_data
   * @return _device_data
   */
  auto& getDeviceData() { return _device_data; }

  /**
   * @brief get _memory_scheduler
   * @return _memory_scheduler
   */
  auto& getMemoryScheduler() const { return _memory_scheduler; }

  void Fill2Device(const std::shared_ptr<MemoryScheduler>& memory_scheduler) {
    this->_memory_scheduler = memory_scheduler;
    this->_memory_scheduler->asyncMemoryH2D();
  }

 private:
  std::shared_ptr<MDData> _md_data;
  std::shared_ptr<PostProcessData> _postProcess_data;
  std::shared_ptr<ConfigData> _config_data;
  std::shared_ptr<DeviceData> _device_data;
  std::shared_ptr<MemoryScheduler> _memory_scheduler;
};
