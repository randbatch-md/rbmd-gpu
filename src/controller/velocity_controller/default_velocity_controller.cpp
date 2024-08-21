#include "default_velocity_controller.h"
#include "velocity_controller_op/update_velocity_op.h"
#include "unit_factor.h"
DefaultVelocityController::DefaultVelocityController() {};

void DefaultVelocityController::Init()
{
    _num_atoms = _structure_info_data->_num_atoms;

    // 配置文件中读取
    _dt = 0.001;           
    auto unit = "LJ";   
    UNIT unit_factor = unit_factor_map[unit]; // 这里可能有重定义隐患

    switch (unit_factor_map["lj"])
    {
    case UNIT::LJ:
        _fmt2v = UnitFactor<UNIT::LJ>::_kb;
    case UNIT::REAL:
        _fmt2v = UnitFactor<UNIT::REAL>::_kb;
    default:
        break;
    }
}

void DefaultVelocityController::Update()
{

    // 这里什么时候用device的position 和 velocity 什么时候用 host的
    bool shake = false;// 配置文件中读取
    if (shake)
    {  
       // 当前LJ不需要添加shake的相关内容
       // __device_data->_shake_vx = _device_data->_d_vx;
       // __device_data->_shake_vy = _device_data->_d_vy;
       // __device_data->_shake_vz = _device_data->_d_vz;
    }

    thrust::raw_point_cat(_device_data->_d_fx.data());

    op::UpdateVelocityOp<device::DEVICE_GPU> UpdateVelocityOp;
    UpdateVelocityOp<device::DEVICE_GPU>(_num_atoms,
                                         _dt, 
                                         _fmt2v, 
                                         _device_data->_d_mass,
                                         _device_data->_d_fx,
                                         _device_data->_d_fy,
                                         _device_data->_d_fz,
                                         _device_data->_d_vx,
                                         _device_data->_d_vy,
                                         _device_data->_d_vz);   
}

