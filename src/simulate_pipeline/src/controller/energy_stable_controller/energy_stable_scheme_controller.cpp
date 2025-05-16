#include "energy_stable_scheme_controller.h"
#include "update_temperature_op.h"

#include <thrust/device_ptr.h>
#include <cmath>

#include "device_types.h"
#include "unit_factor.h"
#include "data_manager.h"

extern rbmd::Real test_e_pe_rbl;
extern rbmd::Real test_e_pe_init;
extern rbmd::Id test_current_step;
EnergyStableSchemeController::EnergyStableSchemeController()  :
  _device_data(DataManager::getInstance().getDeviceData()),
  _structure_info_data(DataManager::getInstance().getMDData()->_structure_info_data)
{
  CHECK_RUNTIME(MALLOC(&_d_temp_contrib, sizeof(rbmd::Real)));
  std::remove("thermo.txt");
  std::remove("temperature.txt");
}
EnergyStableSchemeController::~EnergyStableSchemeController() {
  CHECK_RUNTIME(FREE(_d_temp_contrib));
};

void EnergyStableSchemeController::Init() {
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
    case UNIT::REAL:
      _mvv2e = UnitFactor<UNIT::REAL>::_mvv2e;
      _kB = UnitFactor<UNIT::REAL>::_kb;
      break;
    default:
      break;
  }
}

void EnergyStableSchemeController::Update() {
  ComputeTemp();

  UpdataVelocity();
}

void EnergyStableSchemeController::ComputeTemp() {
  rbmd::Id num_atoms = *(_structure_info_data->_num_atoms);

  CHECK_RUNTIME(MEMSET(_d_temp_contrib, 0, sizeof(rbmd::Real)));

  op::ComputeTemperatureOp<device::DEVICE_GPU>()(
      num_atoms, _mvv2e,
      thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
      thrust::raw_pointer_cast(_device_data->_d_mass.data()),
      thrust::raw_pointer_cast(_device_data->_d_vx.data()),
      thrust::raw_pointer_cast(_device_data->_d_vy.data()),
      thrust::raw_pointer_cast(_device_data->_d_vz.data()), _d_temp_contrib);

  CHECK_RUNTIME(MEMCPY(&_temp_sum, _d_temp_contrib, sizeof(rbmd::Real), D2H));

  bool available_shake = false;

  if (available_shake)
  {
    bool shake = true;
    if (shake) {
      _temperature = 0.5 * _temp_sum / ((3 * num_atoms - num_atoms - 3) * _kB / 2.0);
    } else {
      _temperature = 0.5 * _temp_sum / ((3 * num_atoms - 3) * _kB / 2.0);
    }
  } else
  {
    _temperature = 0.5 * _temp_sum / ((3 * num_atoms - 3) * _kB / 2.0);
  }

  std::cout << "temperature: " << _temperature << std::endl;
  // out
  std::ofstream outfile("temperature.txt", std::ios::app);
  outfile << test_current_step << " " << _temperature << std::endl;
  outfile.close();

  // CHECK_RUNTIME(FREE(temp_contrib));
}

void EnergyStableSchemeController::UpdataVelocity() {
  rbmd::Id num_atoms = *(_structure_info_data->_num_atoms);
  rbmd::Real kinetic_energy = _temp_sum / (2.0*num_atoms);
  rbmd::Real kinetic_energy_init;
  if(1 == test_current_step) {
    kinetic_energy_init = kinetic_energy;
  }

  // Hamiltonian = U_init + kinetic_energy_init;
  rbmd::Real  U_init = test_e_pe_init;
  rbmd::Real  H_init = U_init + kinetic_energy_init;

  //Hamiltonian approximate = U_approximate + kinetic_energy
  rbmd::Real  U_approximate = test_e_pe_rbl;
  rbmd::Real  H_approximate= U_approximate + kinetic_energy ;

  //Rescale  Velocity
  rbmd::Real gamma_ess = 10.0 * _dt;
  rbmd::Real xi_ess =
      SQRT(1.0 +_dt / (gamma_ess*kinetic_energy)* (H_init -H_approximate));
  op::UpdataVelocityRescaleOp<device::DEVICE_GPU> updata_velocity_op;
  updata_velocity_op(num_atoms, xi_ess,
                     thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                     thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                     thrust::raw_pointer_cast(_device_data->_d_vz.data()));

  // out
  std::ofstream outfile("thermo.txt", std::ios::app);
  outfile << test_current_step << " " << _temperature  << " "
     << U_approximate<< " "<<kinetic_energy<< " " <<  H_approximate << std::endl;
  outfile.close();
}
