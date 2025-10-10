#include "default_velocity_controller.h"

#include <thrust/copy.h>
#include <thrust/device_ptr.h>

#include "data_manager.h"
#include "device_types.h"
#include "neighbor_list/include/linked_cell/linked_cell_locator.h"
#include "simulate.h"
#include "unit_factor.h"
#include "update_velocity_op.h"

DefaultVelocityController::DefaultVelocityController(){};

void DefaultVelocityController::Init() {

  auto& num_atoms = *(_structure_info_data->_num_atoms);

  _dt = DataManager::getInstance().getConfigData()->Get<rbmd::Real>(
          "timestep", "execution");//0.001
  auto unit  = DataManager::getInstance().getConfigData()->Get
    <std::string>("unit", "init_configuration", "read_data");
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

  //
  const auto& config = DataManager::getInstance().getConfigData();
  auto integration_type = config->Get<std::string>("integration_type", "execution");
  if("bm" ==integration_type) {
    _par_a = DataManager::getInstance().getConfigData()->Get<rbmd::Real>(
        "par_a", "execution");
    _par_b = DataManager::getInstance().getConfigData()->Get<rbmd::Real>(
            "par_b", "execution");

    _d_prev_fx.resize(num_atoms, 0.0);
    _d_prev_fy.resize(num_atoms, 0.0);
    _d_prev_fz.resize(num_atoms, 0.0);
    _d_pr_prev_fx.resize(num_atoms, 0.0);
    _d_pr_prev_fy.resize(num_atoms, 0.0);
    _d_pr_prev_fz.resize(num_atoms, 0.0);
  }
}

void DefaultVelocityController::Update() {
  bool shake = DataManager::getInstance().getConfigData()->GetJudge<bool>
  ( "fix_shake", "hyper_parameters", "extend");
  if (shake) {
      thrust::copy(_device_data->_d_vx.begin(),_device_data->_d_vx.end(),
        _device_data->_d_shake_vx.begin());
      thrust::copy(_device_data->_d_vy.begin(),_device_data->_d_vy.end(),
        _device_data->_d_shake_vy.begin());
      thrust::copy(_device_data->_d_vz.begin(),_device_data->_d_vz.end(),
        _device_data->_d_shake_vz.begin());
  }

  op::UpdateVelocityOp<device::DEVICE_GPU>()(
      *(_structure_info_data->_num_atoms), _dt, _fmt2v,
      thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
      thrust::raw_pointer_cast(_device_data->_d_mass.data()),
      thrust::raw_pointer_cast(_device_data->_d_fx.data()),
      thrust::raw_pointer_cast(_device_data->_d_fy.data()),
      thrust::raw_pointer_cast(_device_data->_d_fz.data()),
      thrust::raw_pointer_cast(_device_data->_d_vx.data()),
      thrust::raw_pointer_cast(_device_data->_d_vy.data()),
      thrust::raw_pointer_cast(_device_data->_d_vz.data()));

}

void DefaultVelocityController::Updatebm(){
   //std::cout << "test_current_step--v: "  << test_current_step <<  std::endl;
  op::UpdateVelocityOpbm<device::DEVICE_GPU>()(
      *(_structure_info_data->_num_atoms), _par_a,_par_b, _dt, test_current_step, _fmt2v,
      thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
      thrust::raw_pointer_cast(_device_data->_d_mass.data()),
      thrust::raw_pointer_cast(_device_data->_d_fx.data()),
      thrust::raw_pointer_cast(_device_data->_d_fy.data()),
      thrust::raw_pointer_cast(_device_data->_d_fz.data()),
      thrust::raw_pointer_cast(_d_prev_fx.data()),
      thrust::raw_pointer_cast(_d_prev_fy.data()),
      thrust::raw_pointer_cast(_d_prev_fz.data()),
      thrust::raw_pointer_cast(_d_pr_prev_fx.data()),
      thrust::raw_pointer_cast(_d_pr_prev_fy.data()),
      thrust::raw_pointer_cast(_d_pr_prev_fz.data()),
      thrust::raw_pointer_cast(_device_data->_d_vx.data()),
      thrust::raw_pointer_cast(_device_data->_d_vy.data()),
      thrust::raw_pointer_cast(_device_data->_d_vz.data()));
}