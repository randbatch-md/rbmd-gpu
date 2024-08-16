#include "default_velocity_controller.h"
#include "update_velocity_op.h"
DefaultVelocityController::DefaultVelocityController() {};

void DefaultVelocityController::Update()
{
    // 这些都是初始化的参数 看是否需要放到上一层的；
    rbmd::Real dt,_fmt2v;//_dt= 这里需要拿到 time_step _unit_factor._fmt2v

    auto& h_vx = _structure_data->_h_vx;
    auto& h_vy = _structure_data->_h_vy;
    auto& h_vz = _structure_data->_h_vz;

    // 可以直接传入 _force_controller ？
    auto h_fx = _force_controller->_fx;
    auto h_fy = _force_controller->_fy;
    auto h_fz = _force_controller->_fz;

    std::vector<rbmd::Real> mass;
    auto h_mass = _force_field_data->_h_mass; // 这里mass的键值和速度的索引是一一对应的吗

    for (auto it = _h_mass.begin(); it != _h_mass.end(); ++it)
    {
        mass.push_back(it->second);
    }

    _shake_controller->_shake_vx = h_vx;
    _shake_controller->_shake_vy = h_vy;
    _shake_controller->_shake_vz = h_vz;
    _num_atoms
    UpdateVelocityOp<device::DEVICE_GPU>();
    

}

void DefaultVelocityController::Init()
{
    //_dt= 这里需要拿到 time_step
    //_unit_factor._fmt2v




}