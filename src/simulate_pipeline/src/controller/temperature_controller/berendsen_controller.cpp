#include "berendsen_controller.h"

#include <thrust/device_ptr.h>

#include <cmath>

#include "device_types.h"
#include "unit_factor.h"
#include "update_temperature_op.h"
#include "common/thermo_stats.hpp"

rbmd::Real test_temperature;
extern int test_current_step;
BerendsenController::BerendsenController() {
  CHECK_RUNTIME(MALLOC(&_d_temp_contrib, sizeof(rbmd::Real)));
  std::remove("temperature.txt");
}
BerendsenController::~BerendsenController() {
  CHECK_RUNTIME(FREE(_d_temp_contrib));
};

void BerendsenController::Init() {
  auto temperature_array=DataManager::getInstance().getConfigData()->
    GetArray<rbmd::Real>("temperature", "execution"); //[1.0,1.0,0.1]
  _temperature_start = temperature_array[0];
  _temperature_stop = temperature_array[1];
  _temperature_damp = temperature_array[2];

  _dt =  DataManager::getInstance().getConfigData()->Get<rbmd::Real>(
          "timestep", "execution");//0.001

  auto unit = DataManager::getInstance().getConfigData()->Get
    <std::string>("unit", "init_configuration", "read_data");

  UNIT unit_factor = unit_factor_map[unit];
  switch (unit_factor) {
    case UNIT::LJ:
      _mvv2e = UnitFactor<UNIT::LJ>::_mvv2e;
      _kB = UnitFactor<UNIT::LJ>::_kb;
      break;
    case UNIT::METAL:
      _mvv2e = UnitFactor<UNIT::METAL>::_mvv2e;
      _kB = UnitFactor<UNIT::METAL>::_kb;
      break;
    case UNIT::REAL:
      _mvv2e = UnitFactor<UNIT::REAL>::_mvv2e;
      _kB = UnitFactor<UNIT::REAL>::_kb;
      break;
    default:
      break;
  }
  ThermoStats::Instance().AddThermoData("temperature",_temperature_start);
  ThermoStats::Instance().AddThermoData("pressure",0.0);
}

void BerendsenController::Update() {
  //ComputeTemp();

  UpdataVelocity();
}

void BerendsenController::ComputeTemperature() {
  rbmd::Id num_atoms = *(_structure_info_data->_num_atoms);
  CHECK_RUNTIME(MEMSET(_d_temp_contrib, 0, sizeof(rbmd::Real)));

  op::ComputeTemperatureOp<device::DEVICE_GPU>()(num_atoms, _mvv2e,
      thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
      thrust::raw_pointer_cast(_device_data->_d_mass.data()),
      thrust::raw_pointer_cast(_device_data->_d_vx.data()),
      thrust::raw_pointer_cast(_device_data->_d_vy.data()),
      thrust::raw_pointer_cast(_device_data->_d_vz.data()), _d_temp_contrib);

  CHECK_RUNTIME(MEMCPY(&_temp_sum, _d_temp_contrib, sizeof(rbmd::Real), D2H));

  bool available_shake = true; //TODO: need to judge

  if (available_shake)  // H2O
  {
    bool shake = DataManager::getInstance().getConfigData()->GetJudge<bool>
      ("fix_shake", "hyper_parameters", "extend");
    if (shake) {
      _temperature = 0.5 * _temp_sum / ((3 * num_atoms - num_atoms - 3) * _kB / 2.0);
    } else {
      _temperature = 0.5 * _temp_sum / ((3 * num_atoms - 3) * _kB / 2.0);
    }
  } else
  {
    _temperature = 0.5 * _temp_sum / ((3 * num_atoms - 3) * _kB / 2.0);
  }
  test_temperature = _temperature;

  //

  ThermoStats::Instance().AddThermoData("temperature",_temperature);
  ThermoStats::Instance().AddThermoData("pressure",0.0);

  if (std::isnan(_temperature)) {
    Logger::Instance().error( "\033[31mFATAL ERROR: The temperature of the MD simulation is NaN"
                             ". Please check the initial model and the force field parameters. "
    "is invalid.\033[0m");
    exit(EXIT_FAILURE); //
  }

  //out
  auto ensemble_type =DataManager::getInstance().getConfigData()->
    Get<std::string>("ensemble", "execution");

  if ("NVT" == ensemble_type)
  {
    auto interval = DataManager::getInstance().getConfigData()->Get<rbmd::Id>(
"interval", "outputs", "thermo_out");
    std::ofstream outfile("temperature.txt", std::ios::app);
    if (outfile.tellp() == 0) {
      outfile << "step temperature" << std::endl;
    }
    if (test_current_step % interval == 0) {
      outfile << test_current_step << " " << _temperature << std::endl;
    }
    outfile.close();
  }

  // CHECK_RUNTIME(FREE(temp_contrib));
}

void BerendsenController::UpdataVelocity() {
  //
  ComputeTempTargetInit();

  // coeff_berendsen
  rbmd::Real coeff_berendsen =
      SQRT(1.0 + (_dt / _temperature_damp) * (_t_target/ _temperature - 1.0));

  op::UpdataVelocityRescaleOp<device::DEVICE_GPU>()(
                    *(_structure_info_data->_num_atoms), coeff_berendsen,
                     thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                     thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                     thrust::raw_pointer_cast(_device_data->_d_vz.data()));
}
