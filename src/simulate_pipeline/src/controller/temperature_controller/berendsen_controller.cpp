#include "berendsen_controller.h"
#include "update_temperature_op.h"
#include "unit_factor.h"
#include "device_types.h"
#include <thrust/device_ptr.h>
#include <cmath>

BerendsenController::BerendsenController() {};

void BerendsenController::Init()
{
    _num_atoms = *(_structure_info_data->_num_atoms);
    _temp_sum = 0;
    _Tdamp = 0.1;    // 配置文件中读取 temperature [1.0,1.0,0.1]
    _dt = 0.001;     // 配置文件中读取      
    auto unit = "LJ";   
    UNIT unit_factor = unit_factor_map[unit]; // 这里可能有重定义隐患

    switch (unit_factor)
    {
    case UNIT::LJ:
        _mvv2e = UnitFactor<UNIT::LJ>::_mvv2e;
        _kB = UnitFactor<UNIT::LJ>::_kb;
        break;
    case UNIT::REAL:
        _mvv2e = UnitFactor<UNIT::REAL>::_mvv2e;
        _kB = UnitFactor<UNIT::REAL>::_kb;
        break;
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
    //rbmd::Real* temp_contrib;
    //CHECK_RUNTIME(hipMemset(temp_contrib, 0, sizeof(rbmd::Real)));

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

    bool available_shake = false;

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

void BerendsenController::UpdataVelocity()
{
    rbmd::Real kbT = 1;//配置文件获取 
    rbmd::Real coeff_Berendsen = std::sqrt(1.0 + (_dt / _Tdamp) * (kbT / _temp - 1.0)); 

    op::UpdataVelocityRescaleOp<device::DEVICE_GPU> updata_velocity_op; 
    updata_velocity_op(_num_atoms,
                       coeff_Berendsen,
                       thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                       thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                       thrust::raw_pointer_cast(_device_data->_d_vz.data()));
}
