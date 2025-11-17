#pragma once

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include "data_manager.h"
#include "model/device_data.h"
#include "model/md_data.h"
#include "model/structure_info_data.h"

class GroupController {
public:

  GroupController();
  ~GroupController();
  GroupController(const GroupController&) = delete;
  GroupController& operator=(const GroupController&) = delete;

  // 初始化，通常在模拟开始时调用一次
  void Init();

  /**
   * @brief 计算指定原子组的质心速度 (Velocity Center-of-Mass)
   * @param group_name 要计算的组名
   * @param vcm_out 用于存储计算结果的数组 (大小为3的host端数组)
   */
  void ComputeVCM(const std::string& group_name, Real3 vcm_out);


 /**
  * @brief 计算指定原子组的质心位置 (Center-of-Mass Position)
  * @param group_name 要计算的组名
  * @param xcm_out 用于存储计算结果的数组 (大小为3的host端数组)
  */
 void ComputeXCM(const std::string& group_name, Real3 xcm_out);

    /**
     * @brief 计算指定原子组相对于某点的角动量
     * @param group_name 要计算的组名
     * @param origin 计算角动量的参考原点 (通常是质心xcm)
     * @param angmom_out 用于存储计算结果的数组 (大小为3的host端数组)
     */
    void ComputeAngMom(const std::string& group_name, Real3 cm, rbmd::Real* angmom_out);

    /**
     * @brief 计算指定原子组相对于某点的转动惯量张量
     * @param group_name 要计算的组名
     * @param cm 参考原点 (通常是质心xcm)
     * @param inertia_out 用于存储计算结果的二维数组 (大小为3x3的host端数组)
     */
    void ComputeInertia(const std::string& group_name, Real3 origin, rbmd::Real (*inertia_out)[3]);

    /**
     * @brief 根据角动量和转动惯量计算角速度 (在CPU上完成),
     * @param angmom 角动量 L
     * @param inertia 转动惯量张量 I
     * @param omega_out 计算出的角速度 omega = I^-1 * L
     */
    void ComputeOmega(const rbmd::Real* angmom, const rbmd::Real (*inertia)[3], rbmd::Real* omega_out);

private:

  // 缓存数据指针，避免反复调用DataManager
  std::unordered_map<std::string, std::vector<rbmd::Id>> _groups;
  std::shared_ptr<DeviceData> _device_data;
  std::shared_ptr<StructureInfoData> _structure_info_data;
  std::shared_ptr<Box> _box;

};