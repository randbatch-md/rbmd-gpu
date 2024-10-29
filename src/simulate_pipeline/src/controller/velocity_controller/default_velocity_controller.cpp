#include "default_velocity_controller.h"

#include <thrust/device_ptr.h>

#include "data_manager.h"
#include "device_types.h"
#include "neighbor_list/include/linked_cell/linked_cell_locator.h"
#include "unit_factor.h"
#include "update_velocity_op.h"
#include <thrust/copy.h>

DefaultVelocityController::DefaultVelocityController(){};

void DefaultVelocityController::Init() {

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
}

void DefaultVelocityController::Update() {
  bool shake = DataManager::getInstance().getConfigData()->GetJudge<bool>( "fix_shake", "hyper_parameters", "extend");
  if (shake) {
      thrust::copy(_device_data->_d_vx.begin(), _device_data->_d_vx.end(), _device_data->_d_shake_vx.begin());
      thrust::copy(_device_data->_d_vy.begin(), _device_data->_d_vy.end(), _device_data->_d_shake_vy.begin());
      thrust::copy(_device_data->_d_vz.begin(), _device_data->_d_vz.end(), _device_data->_d_shake_vz.begin());

      std::cout<<"---------更新速度内部--------------"<<std::endl;
      std::vector<rbmd::Real> h_vx(300);
      std::vector<rbmd::Real> h_px(300);
      std::vector<rbmd::Real> h_shakeA_vx(300);
      std::vector<rbmd::Real> h_shakeA_px(300);
      thrust::copy(_device_data->_d_vx.begin(), _device_data->_d_vx.end(), h_vx.begin());
      thrust::copy(_device_data->_d_px.begin(), _device_data->_d_px.end(), h_px.begin());
      thrust::copy(_device_data->_d_shake_vx.begin(), _device_data->_d_shake_vx.end(), h_shakeA_vx.begin());
      thrust::copy(_device_data->_d_shake_px.begin(), _device_data->_d_shake_px.end(), h_shakeA_px.begin());

      for (int j = 0; j < 10; ++j) {std::cout<<h_vx[j]<<" , ";}
      std::cout<<std::endl;
      for (int j = 0; j < 10; ++j) {std::cout<<h_px[j]<<" , ";}
      std::cout<<std::endl;

      for (int j = 0; j < 10; ++j) {std::cout<< h_shakeA_vx[j]<<" , ";}
      std::cout<<std::endl;
      for (int j = 0; j < 10; ++j) {std::cout<< h_shakeA_px[j]<<" , ";}
      std::cout<<std::endl;
  }

  op::UpdateVelocityOp<device::DEVICE_GPU> update_velocity_op;
  update_velocity_op(
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
