#include "default_velocity_controller.h"

#include <thrust/device_ptr.h>
#include <thrust/copy.h>
#include "data_manager.h"
#include "device_types.h"
#include "neighbor_list/include/linked_cell/linked_cell_locator.h"
#include "unit_factor.h"
#include "update_velocity_op.h"

DefaultVelocityController::DefaultVelocityController(){};

void DefaultVelocityController::Init() {
  _num_atoms = *(_structure_info_data->_num_atoms);

  _dt = DataManager::getInstance().getConfigData()->Get<rbmd::Real>( "timestep", "execution"); // TODO: Json file
  auto unit = DataManager::getInstance().getConfigData()->Get<std::string>( "unit", "init_configuration", "read_data");

  UNIT unit_factor = unit_factor_map[unit];

  switch (unit_factor) {
    case UNIT::METAL:
      _fmt2v = UnitFactor<UNIT::METAL>::_fmt2v;
      break;
    case UNIT::LJ:
      _fmt2v = UnitFactor<UNIT::LJ>::_fmt2v;
      break;
    case UNIT::REAL:
      _fmt2v = UnitFactor<UNIT::REAL>::_fmt2v;
      break;
    default:
      break;
  }
}

void DefaultVelocityController::Update() {
  bool shake = DataManager::getInstance().getConfigData()->GetJudge<bool>( "fix_shake", "hyper_parameters", "extend");
  if (shake) {
      thrust::copy(_device_data->_d_vx.begin(), _device_data->_d_vx.end(), _device_data->_d_shake_vx.begin());
      thrust::copy(_device_data->_d_vy.begin(), _device_data->_d_vy.end(), _device_data->_d_shake_vy.begin());
      thrust::copy(_device_data->_d_vz.begin(), _device_data->_d_vz.end(), _device_data->_d_shake_vz.begin());
      //thrust::host_vector<rbmd::Real> host_shake_vx = _device_data->_d_shake_vx;
      //for (const auto& value : host_shake_vx) {std::cout << value << " ";}
  }

  op::UpdateVelocityOp<device::DEVICE_GPU> update_velocity_op;
  update_velocity_op(
      _num_atoms, _dt, _fmt2v,
      thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
      thrust::raw_pointer_cast(_device_data->_d_mass.data()),
      thrust::raw_pointer_cast(_device_data->_d_fx.data()),
      thrust::raw_pointer_cast(_device_data->_d_fy.data()),
      thrust::raw_pointer_cast(_device_data->_d_fz.data()),
      thrust::raw_pointer_cast(_device_data->_d_vx.data()),
      thrust::raw_pointer_cast(_device_data->_d_vy.data()),
      thrust::raw_pointer_cast(_device_data->_d_vz.data()));
}
