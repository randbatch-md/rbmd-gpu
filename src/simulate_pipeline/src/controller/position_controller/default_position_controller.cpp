#include "default_position_controller.h"

#include <thrust/device_ptr.h>

#include "neighbor_list/include/linked_cell/linked_cell_locator.h"
#include "update_position_op.h"

DefaultPositionController::DefaultPositionController(){};

void DefaultPositionController::Init() 
{
  _dt = DataManager::getInstance().getConfigData()->Get<rbmd::Real>("timestep", "execution");//0.001
}

void DefaultPositionController::Update() {
  bool available_shake = false;
  if (available_shake) {
    bool shake = false;
    if (shake) {
      //_device_data->_shake_vx = _device_data->_d_px;
      //_device_data->_shake_vy = _device_data->_d_py;
      //_device_data->_shake_vz = _device_data->_d_pz;
    }
    op::UpdatePositionOp<device::DEVICE_GPU>()(*(_structure_info_data->_num_atoms), _dt, _device_data->_d_box,
                       thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                       thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                       thrust::raw_pointer_cast(_device_data->_d_vz.data()),
                       thrust::raw_pointer_cast(_device_data->_d_px.data()),
                       thrust::raw_pointer_cast(_device_data->_d_py.data()),
                       thrust::raw_pointer_cast(_device_data->_d_pz.data()));
  } else {
    op::UpdatePositionFlagOp<device::DEVICE_GPU>()(*(_structure_info_data->_num_atoms), _dt, _device_data->_d_box,
                       thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                       thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                       thrust::raw_pointer_cast(_device_data->_d_vz.data()),
                       thrust::raw_pointer_cast(_device_data->_d_px.data()),
                       thrust::raw_pointer_cast(_device_data->_d_py.data()),
                       thrust::raw_pointer_cast(_device_data->_d_pz.data()),
                       thrust::raw_pointer_cast(_device_data->_d_flagX.data()),
                       thrust::raw_pointer_cast(_device_data->_d_flagY.data()),
                       thrust::raw_pointer_cast(_device_data->_d_flagZ.data()));
  }
}

void DefaultPositionController::SetCenterTargetPositions() {
  std::string init_type = "inbuild";
  if (init_type == _init_type) {
  }
}
