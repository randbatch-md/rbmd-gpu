#include "default_velocity_controller.h"
#include "velocity_controller_op/update_velocity_op.h"
#include "unit_factor.h"
DefaultVelocityController::DefaultVelocityController() {};

void DefaultVelocityController::Update()
{

    // ����ʲôʱ����device��position �� velocity ʲôʱ���� host��
    bool shake = false;// �����ļ��ж�ȡ
    if (shake)
    {
        // ��Щ���ǳ�ʼ���Ĳ��� ���Ƿ���Ҫ�ŵ���һ��ģ�
        auto& h_vx = _structure_data->_h_vx;
        auto& h_vy = _structure_data->_h_vy;
        auto& h_vz = _structure_data->_h_vz;

        // ����ֱ�Ӵ��� _force_controller ��    
        _shake_controller->_shake_vx = h_vx;
        _shake_controller->_shake_vy = h_vy;
        _shake_controller->_shake_vz = h_vz;
    }

    op::UpdateVelocityOp<device::DEVICE_GPU> UpdateVelocityOp;
    UpdateVelocityOp<device::DEVICE_GPU>(_structure_info_data->_num_atoms,
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

void DefaultVelocityController::Init()
{
    _dt=0.001           // �����ļ��ж�ȡ
    auto unit = "LJ";   // �����ļ��ж�ȡ
    UNIT unit_factor = unit_factor_map[unit]; // ����������ض�������

    switch (unit_factor_map["lj"])
    {
    case UNIT::LJ:
        _fmt2v = UnitFactor<UNIT::LJ>::_kb;
    case UNIT::REAL:
        _fmt2v = UnitFactor<UNIT::REAL>::_kb;
    default:
        break;
    }
    _num_atoms = _structure_info_data->_num_atoms;
}