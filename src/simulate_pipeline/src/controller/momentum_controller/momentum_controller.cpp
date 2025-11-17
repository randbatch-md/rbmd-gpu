#include "momentum_controller.h"

#include "../group_controller/group_controller.h"
#include "common/thermo_stats.hpp"
#include "momentum_controller_op.h"
#include "unit_factor.h"
#include "update_temperature_op.h"
extern int test_current_step;

MomentumController::MomentumController()
    : _interval(1), _linear_flag(false), _angular_flag(false), _rescale_flag(false),
      _x_flag(false), _y_flag(false), _z_flag(false), _group_name("all")
{
    _structure_info_data = DataManager::getInstance().getMDData()->_structure_info_data;
    _device_data = DataManager::getInstance().getDeviceData();
    _box = DataManager::getInstance().getMDData()->_box;

    _group_controller = std::make_shared<GroupController>();
    CHECK_RUNTIME(MALLOC(&_d_ke_contrib, sizeof(rbmd::Real)));

  auto unit = DataManager::getInstance().getConfigData()->Get
  <std::string>("unit", "init_configuration", "read_data");

  UNIT unit_factor = unit_factor_map[unit];
  switch (unit_factor) {
    case UNIT::LJ:
      _mvv2e = UnitFactor<UNIT::LJ>::_mvv2e;
      break;
    case UNIT::METAL:
      _mvv2e = UnitFactor<UNIT::METAL>::_mvv2e;
      break;
    case UNIT::REAL:
      _mvv2e = UnitFactor<UNIT::REAL>::_mvv2e;
      break;
    default:
      break;
  }
}

MomentumController::~MomentumController() {
    CHECK_RUNTIME(FREE(_d_ke_contrib));
}

void MomentumController::Init() {
    auto config = DataManager::getInstance().getConfigData();

    // 检查是否启用了动量控制
    if (!config->PathExists({"execution", "momentum_control"}))
      return;

    _interval = config->Get<int>("interval", "execution", "momentum_control");
    
    // 读取linear关键字及其标志
  if (config->PathExists({"execution", "momentum_control","linear"})) {
        _linear_flag = true;
        auto linear_flags = config->GetArray<int>("linear", "execution", "momentum_control");
        if (linear_flags.size() == 3) {
            _x_flag = static_cast<bool>(linear_flags[0]);
            _y_flag = static_cast<bool>(linear_flags[1]);
            _z_flag = static_cast<bool>(linear_flags[2]);
        }
    }

    // 读取angular关键字
    _angular_flag = config->GetJudge<bool>("angular", "execution", "momentum_control");

    // 读取rescale关键字
    _rescale_flag = config->GetJudge<bool>("rescale", "execution", "momentum_control");
    
    _group_name = config->Get<std::string>("group", "execution", "momentum_control");
    if (_group_name.empty()) {
        _group_name = "all";
    }

}

void MomentumController::Execute() {
    if (_interval == 0) return;

    if (test_current_step % _interval != 0)
      return;

    rbmd::Real ekin_old = 1.0;
    rbmd::Real ekin_new = 1.0;
    rbmd::Id num_atoms = *(_structure_info_data->_num_atoms);

    // 步骤1: (可选) 计算并存储原始动能
    if (_rescale_flag)
    {
        thrust::device_vector<rbmd::Real> d_ke_contrib(1, 1.0);
        op::ComputeTemperatureOp<device::DEVICE_GPU>()(num_atoms, _mvv2e,
      thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
      thrust::raw_pointer_cast(_device_data->_d_mass.data()),
      thrust::raw_pointer_cast(_device_data->_d_vx.data()),
      thrust::raw_pointer_cast(_device_data->_d_vy.data()),
      thrust::raw_pointer_cast(_device_data->_d_vz.data()),
      thrust::raw_pointer_cast(d_ke_contrib.data()));

       CHECK_RUNTIME(MEMCPY(&ekin_old, thrust::raw_pointer_cast(d_ke_contrib.data()),
          sizeof(rbmd::Real), D2H));
    }
    
    // 步骤2: 移除线动量
    if (_linear_flag)
    {
        Real3 vcm{0.0,0.0,0.0};
       _group_controller->ComputeVCM(_group_name, vcm);

        op::ZeroLinearMomentumOp<device::DEVICE_GPU>()(
            num_atoms, vcm, _x_flag, _y_flag, _z_flag,
            thrust::raw_pointer_cast(_device_data->_d_vx.data()),
            thrust::raw_pointer_cast(_device_data->_d_vy.data()),
            thrust::raw_pointer_cast(_device_data->_d_vz.data())
        );

    }

    // 步骤3: 移除角动量
    if (_angular_flag)
    {
        Real3 xcm{0.0,0.0,0.0};
        rbmd::Real angmom[3], omega[3];
        rbmd::Real inertia[3][3];
        
        _group_controller->ComputeXCM(_group_name, xcm);
        _group_controller->ComputeAngMom(_group_name, xcm, angmom);
        _group_controller->ComputeInertia(_group_name, xcm, inertia);

        // 计算角速度 omega = I^-1 * L
        // 这个3x3矩阵求逆和乘法在CPU上完成即可
        _group_controller->ComputeOmega(angmom, inertia, omega);

      thrust::device_vector<rbmd::Real> d_omega(3, 0.0);
      CHECK_RUNTIME(MEMCPY(thrust::raw_pointer_cast(d_omega.data()), omega,
        3*sizeof(rbmd::Real), H2D));

      op::ZeroAngularMomentumOp<device::DEVICE_GPU>()(
            num_atoms, xcm, thrust::raw_pointer_cast(d_omega.data()),
            thrust::raw_pointer_cast(_device_data->_d_px.data()),
            thrust::raw_pointer_cast(_device_data->_d_py.data()),
            thrust::raw_pointer_cast(_device_data->_d_pz.data()),
            thrust::raw_pointer_cast(_device_data->_d_flagX.data()),
            thrust::raw_pointer_cast(_device_data->_d_flagY.data()),
            thrust::raw_pointer_cast(_device_data->_d_flagZ.data()),
            *_box,
            thrust::raw_pointer_cast(_device_data->_d_vx.data()),
            thrust::raw_pointer_cast(_device_data->_d_vy.data()),
            thrust::raw_pointer_cast(_device_data->_d_vz.data()));
    }

    // 步骤4: (可选) 重新缩放速度以恢复动能
    if (_rescale_flag)
    {
      thrust::device_vector<rbmd::Real> d_ke_contrib_2(1, 1.0);
        op::ComputeTemperatureOp<device::DEVICE_GPU>()(num_atoms, _mvv2e,
  thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
  thrust::raw_pointer_cast(_device_data->_d_mass.data()),
  thrust::raw_pointer_cast(_device_data->_d_vx.data()),
  thrust::raw_pointer_cast(_device_data->_d_vy.data()),
  thrust::raw_pointer_cast(_device_data->_d_vz.data()),
  thrust::raw_pointer_cast(d_ke_contrib_2.data()));

      //
      CHECK_RUNTIME(MEMCPY(&ekin_new, thrust::raw_pointer_cast(d_ke_contrib_2.data()),
        sizeof(rbmd::Real), D2H));
      rbmd::Real factor = 1.0;
      if (ekin_new > 1e-9)
        factor = SQRT(ekin_old / ekin_new);

      // std::cout <<  "factor: " << factor <<", ekin_new: "
      //   <<  ekin_new << ", ekin_new: "<< ekin_new<<  std::endl;

      op::UpdataVelocityRescaleOp<device::DEVICE_GPU>()(
            num_atoms, factor,
            thrust::raw_pointer_cast(_device_data->_d_vx.data()),
            thrust::raw_pointer_cast(_device_data->_d_vy.data()),
            thrust::raw_pointer_cast(_device_data->_d_vz.data()));
  }

}