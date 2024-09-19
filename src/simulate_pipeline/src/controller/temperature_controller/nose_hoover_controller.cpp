#include "nose_hoover_controller.h"
#include "update_temperature_op.h"
#include "unit_factor.h"
#include "device_types.h"
#include <thrust/device_ptr.h>
#include <cmath>

NoseHooverController::NoseHooverController() {};

void NoseHooverController::Init()
{
    _num_atoms = *(_structure_info_data->_num_atoms);
    _temp_sum = 0;
    _nosehooverxi = 0;
    _dt = 0.001;     // 配置文件中读取      
    auto unit = "LJ";   
    UNIT unit_factor = unit_factor_map[unit]; // 这里可能有重定义隐患

    switch (unit_factor)
    {
    case UNIT::LJ:
        _mvv2e = UnitFactor<UNIT::LJ>::_mvv2e;
        _kB = UnitFactor<UNIT::LJ>::_kb;
        _fmt2v = UnitFactor<UNIT::LJ>::_fmt2v;
        break;
    case UNIT::REAL:
        _mvv2e = UnitFactor<UNIT::REAL>::_mvv2e;
        _kB = UnitFactor<UNIT::REAL>::_kb;
        _fmt2v = UnitFactor<UNIT::REAL>::_fmt2v;
        break;

    default:
        break;
    }
}

void NoseHooverController::Update()
{
    ComputeTemp();

    UpdataVelocity();
}

void NoseHooverController::ComputeTemp()
{
    _temp_sum = 0;

    rbmd::Real* temp_contrib;
    CHECK_RUNTIME(MALLOC(&temp_contrib, sizeof(rbmd::Real)));
    CHECK_RUNTIME(MEMCPY(temp_contrib, &_temp_sum, sizeof(rbmd::Real), H2D));

    op::ComputeTemperatureOp<device::DEVICE_GPU> compute_temperature_op;
    compute_temperature_op(_num_atoms,
                           _mvv2e,
                           thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
                           thrust::raw_pointer_cast(_device_data->_d_mass.data()),
                           thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                           thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                           thrust::raw_pointer_cast(_device_data->_d_vz.data()),
                           temp_contrib);

    CHECK_RUNTIME(MEMCPY(&_temp_sum, temp_contrib, sizeof(rbmd::Real), D2H));

    bool available_shake = true;

    if (available_shake) // H2O / NACl / EAM ...
    {
        bool shake = true;
        if (shake)
        {
            _temp = 0.5 * _temp_sum / ((3 * _num_atoms - _num_atoms - 3) * _kB / 2.0);
        }
        else
        {
            _temp = 0.5 * _temp_sum / ((3 * _num_atoms - 3) * _kB / 2.0);
        }
    }
    else // PEO 
    {
        _temp = 0.5 * _temp_sum / ((3 * _num_atoms - 3) * _kB / 2.0);

    }

    std::cout << "_temp=" << _temp << std::endl;
    CHECK_RUNTIME(FREE(temp_contrib));
}

void NoseHooverController::UpdataVelocity()
{
    rbmd::Real kbT = 1;//配置文件获取 

    op::UpdataVelocityNoseHooverOp<device::DEVICE_GPU> updata_velocity_nose_hoover_op;
    updata_velocity_nose_hoover_op(_num_atoms,
                                   _dt,
                                   _fmt2v,
                                   _nosehooverxi,
                                   thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
                                   thrust::raw_pointer_cast(_device_data->_d_mass.data()),
                                   thrust::raw_pointer_cast(_device_data->_d_fx.data()),
                                   thrust::raw_pointer_cast(_device_data->_d_fy.data()),
                                   thrust::raw_pointer_cast(_device_data->_d_fz.data()),
                                   thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                                   thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                                   thrust::raw_pointer_cast(_device_data->_d_vz.data()));


    _nosehooverxi += 0.5 * _dt * (_temp / kbT - 1.0) / (std::pow(10.0, -1) * _dt);

}
