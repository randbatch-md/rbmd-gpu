#include "berendsen_controller.h"
#include "update_temperature_op.h"
#include "unit_factor.h"
#include "device_types.h"
#include <thrust/device_ptr.h>
#include <cmath>

BerendsenController::BerendsenController() {};

void BerendsenController::Init()
{
    _num_atoms = _structure_info_data->_num_atoms;
    _temp_sum = 0;
    _Tdamp = 0.1;    // �����ļ��ж�ȡ temperature [1.0,1.0,0.1]
    _dt = 0.001;     // �����ļ��ж�ȡ      
    auto unit = "LJ";   
    UNIT unit_factor = unit_factor_map[unit]; // ����������ض�������

    switch (unit_factor)
    {
    case UNIT::LJ:
        _mvv2e = UnitFactor<UNIT::LJ>::_mvv2e;
        _kB = UnitFactor<UNIT::LJ>::_kb;
    case UNIT::REAL:
        _mvv2e = UnitFactor<UNIT::REAL>::_mvv2e;
        _kB = UnitFactor<UNIT::REAL>::_kb;
    default:
        break;
    }
}

void BerendsenController::Update()
{
    ComputeTemp();

    UpdataVelocity();
}

void BerendsenController::ComputeTemp()
{
    op::ComputeTemperatureOp<device::DEVICE_GPU> compute_temperature_op;
    compute_temperature_op(_num_atoms,
                           _mvv2e,
                           thrust::raw_pointer_cast(_device_data->_d_mass.data()),
                           thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                           thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                           thrust::raw_pointer_cast(_device_data->_d_vz.data()),
                           _temp_sum);

    bool available_shake = false;
    std::string init_type = "inbuild";
    if (init_type== "inbuild")  // LJ / LJ_salt
    {
        _temp = 0.5 * _temp_sum / (3 * _num_atoms / 2.0);
    }
    else
    {
        if (available_shake) // H2O / NACl / EAM ...
        {
            bool shake = false;
            if (shake)
            {
                _temp = 0.5 * _temp_sum / ((3 * _num_atoms - _num_atoms - 3) * _kB / 2.0);
            }
            else
            {
                _temp = 0.5 * _temp_sum / ((3 * _num_atoms) * _kB / 2.0);
            }
        }
        else // PEO
        {
            _temp = 0.5 * _temp_sum / ((3 * _num_atoms - 3) * _kB / 2.0);

        }
    }
}

void BerendsenController::UpdataVelocity()
{
    rbmd::Real kbT = 1;//�����ļ���ȡ 
    rbmd::Real coeff_Berendsen = std::sqrt(1.0 + (_dt / _Tdamp) * (kbT / _temp - 1.0)); 

    op::UpdataVelocityRescaleOp<device::DEVICE_GPU> updata_velocity_op; 
    updata_velocity_op(_num_atoms,
                       coeff_Berendsen,
                       thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                       thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                       thrust::raw_pointer_cast(_device_data->_d_vz.data()));
}