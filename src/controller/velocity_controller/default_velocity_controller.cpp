#include "default_velocity_controller.h"
#include "update_velocity_op.h"
DefaultVelocityController::DefaultVelocityController() {};

void DefaultVelocityController::Update()
{
    // ��Щ���ǳ�ʼ���Ĳ��� ���Ƿ���Ҫ�ŵ���һ��ģ�
    rbmd::Real dt,_fmt2v;//_dt= ������Ҫ�õ� time_step _unit_factor._fmt2v

    auto& h_vx = _structure_data->_h_vx;
    auto& h_vy = _structure_data->_h_vy;
    auto& h_vz = _structure_data->_h_vz;

    // ����ֱ�Ӵ��� _force_controller ��
    auto h_fx = _force_controller->_fx;
    auto h_fy = _force_controller->_fy;
    auto h_fz = _force_controller->_fz;

    std::vector<rbmd::Real> mass;
    auto h_mass = _force_field_data->_h_mass; // ����mass�ļ�ֵ���ٶȵ�������һһ��Ӧ����

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
    //_dt= ������Ҫ�õ� time_step
    //_unit_factor._fmt2v




}