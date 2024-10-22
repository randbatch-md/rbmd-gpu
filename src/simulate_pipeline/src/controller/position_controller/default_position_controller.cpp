#include "default_position_controller.h"

#include <thrust/device_ptr.h>

#include "neighbor_list/include/linked_cell/linked_cell_locator.h"
#include "update_position_op.h"

DefaultPositionController::DefaultPositionController(){};

void DefaultPositionController::Init() {
  // auto atom_style = _config_data->Get<std::string>("atom_style",
  // "init_configuration", "read_data");

  _dt = DataManager::getInstance().getConfigData()->Get<rbmd::Real>( "timestep", "execution"); // TODO: Json file
  _num_atoms = *(_structure_info_data->_num_atoms);
}

void DefaultPositionController::Update() {
  bool available_shake = false;//TODO: json file to judge
  if (available_shake) {
    bool shake = DataManager::getInstance().getConfigData()->GetJudge<bool>( "fix_shake", "hyper_parameters", "extend");
    if (shake) {
        thrust::copy(_device_data->_d_px.begin(), _device_data->_d_px.end(), _device_data->_d_px.begin());
        thrust::copy(_device_data->_d_py.begin(), _device_data->_d_py.end(), _device_data->_d_py.begin());
        thrust::copy(_device_data->_d_pz.begin(), _device_data->_d_pz.end(), _device_data->_d_pz.begin());
    }
    op::UpdatePositionOp<device::DEVICE_GPU> update_position_op;
    update_position_op(_num_atoms, _dt, (*_structure_info_data->_range)[0][0],
                       (*_structure_info_data->_range)[1][0],
                       (*_structure_info_data->_range)[2][0],
                       (*_structure_info_data->_range)[0][1],
                       (*_structure_info_data->_range)[1][1],
                       (*_structure_info_data->_range)[2][1],
                       thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                       thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                       thrust::raw_pointer_cast(_device_data->_d_vz.data()),
                       thrust::raw_pointer_cast(_device_data->_d_px.data()),
                       thrust::raw_pointer_cast(_device_data->_d_py.data()),
                       thrust::raw_pointer_cast(_device_data->_d_pz.data()));
  }
    else
    {
    op::UpdatePositionFlagOp<device::DEVICE_GPU> update_position_op;
    update_position_op(_num_atoms, _dt, _device_data->_d_box,
                       (*_structure_info_data->_range)[0][0],
                       (*_structure_info_data->_range)[1][0],
                       (*_structure_info_data->_range)[2][0],
                       (*_structure_info_data->_range)[0][1],
                       (*_structure_info_data->_range)[1][1],
                       (*_structure_info_data->_range)[2][1],
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
