#include "rescale_controller.h"
#include "temperature_controller_op/update_temperature_op.h"
#include "../common/unit_factor.h"
#include <thrust/device_ptr.h>

RescaleController::RescaleController() {};

void RescaleController::Init()
{
    _num_atoms = _structure_info_data->_num_atoms;
    _temp_sum = 0;
      
    auto unit = "LJ";   //配置文件读取
    UNIT unit_factor = unit_factor_map[unit]; // 这里可能有重定义隐患

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

void RescaleController::Update()
{
    ComputeTemp();

    UpdataVelocity();
}

void RescaleController::ComputeTemp()
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
    if (init_type== "inbuild") // LJ / LJ_salt
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

void RescaleController::UpdataVelocity()
{
    rbmd::Real kbT = 1;//配置文件获取
    rbmd::Real coeff_rescale = std::sqrt(kbT / _temp);

    
    op::UpdataVelocityRescaleOp<device::DEVICE_GPU> updata_velocity_op;
    updata_velocity_op(_num_atoms,
                       coeff_rescale,
                       thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                       thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                       thrust::raw_pointer_cast(_device_data->_d_vz.data()));
}
