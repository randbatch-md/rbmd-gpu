#pragma once

#include <memory>
#include <string>

#include "data_manager.h"
#include "group_controller.h"
#include "model/md_data.h"

class MomentumController {
public:
  MomentumController();
  virtual ~MomentumController();

  /**
   * @brief 从配置文件读取参数并初始化
   */
  void Init();

  /**
   * @brief 在每个时间步末尾执行动量修正
   */
  void Execute();

private:
  std::shared_ptr<DeviceData> _device_data;
  std::shared_ptr<StructureInfoData> _structure_info_data;
  std::shared_ptr<Box> _box;

  // --- 从配置文件读取的参数 ---
  int _interval;
  bool _linear_flag;
  bool _angular_flag;
  bool _rescale_flag;
  bool _x_flag, _y_flag, _z_flag;
  std::string _group_name;

  // --- GPU端用于中间计算的内存 ---
  rbmd::Real* _d_ke_contrib; // 用于计算动能


  std::shared_ptr<GroupController> _group_controller;

  rbmd::Real _mvv2e;
};