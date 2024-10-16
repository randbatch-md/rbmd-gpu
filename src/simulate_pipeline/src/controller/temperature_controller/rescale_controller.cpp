#include "rescale_controller.h"

#include <thrust/device_ptr.h>

#include "device_types.h"
#include "neighbor_list/include/linked_cell/linked_cell_locator.h"
#include "rbmd_define.h"
#include "unit_factor.h"
#include "update_temperature_op.h"

RescaleController::RescaleController() {
  CHECK_RUNTIME(MALLOC(&_d_temp_contrib, sizeof(rbmd::Real)));
}
RescaleController::~RescaleController() {
  CHECK_RUNTIME(FREE(_d_temp_contrib));
}

void RescaleController::Init() {
  _num_atoms = *(_structure_info_data->_num_atoms);
  _temp_sum = 0;

  auto unit = "LJ";
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

void RescaleController::Update() {
  ComputeTemp();

  UpdataVelocity();
}

void RescaleController::ComputeTemp() {
  extern int test_current_step;
  // CHECK_RUNTIME(MEMCPY(temp_contrib, &_temp_sum, sizeof(rbmd::Real), H2D));
  CHECK_RUNTIME(MEMSET(_d_temp_contrib, 0, sizeof(rbmd::Real)));
  op::ComputeTemperatureOp<device::DEVICE_GPU> compute_temperature_op;
  compute_temperature_op(
      _num_atoms, _mvv2e,
      thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
      thrust::raw_pointer_cast(_device_data->_d_mass.data()),
      thrust::raw_pointer_cast(_device_data->_d_vx.data()),
      thrust::raw_pointer_cast(_device_data->_d_vy.data()),
      thrust::raw_pointer_cast(_device_data->_d_vz.data()), _d_temp_contrib);

  CHECK_RUNTIME(MEMCPY(&_temp_sum, _d_temp_contrib, sizeof(rbmd::Real), D2H));
  // std::cout << "_temp_sum" << _temp_sum << std::endl;

  bool available_shake = false;

  if (available_shake)  // H2O / NACl / EAM ...
  {
    bool shake = true;
    if (shake) {
      _temp = 0.5 * _temp_sum / ((3 * _num_atoms - _num_atoms - 3) * _kB / 2.0);
    } else {
      _temp = 0.5 * _temp_sum / ((3 * _num_atoms - 3) * _kB / 2.0);
    }
  } else  // PEO
  {
    _temp = 0.5 * _temp_sum / ((3 * _num_atoms - 3) * _kB / 2.0);
  }

  std::cout <<  "step=" << test_current_step <<
    ", " << "_temp=" << _temp << std::endl;

  // out
  std::ofstream outfile("temp.txt", std::ios::app);
  outfile << test_current_step << " " << _temp << std::endl;
  outfile.close();
}

void RescaleController::UpdataVelocity() {
  rbmd::Real kbT = 1;
  rbmd::Real coeff_rescale = std::sqrt(kbT / _temp);

  op::UpdataVelocityRescaleOp<device::DEVICE_GPU> updata_velocity_op;
  updata_velocity_op(_num_atoms, coeff_rescale,
                     thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                     thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                     thrust::raw_pointer_cast(_device_data->_d_vz.data()));
}
