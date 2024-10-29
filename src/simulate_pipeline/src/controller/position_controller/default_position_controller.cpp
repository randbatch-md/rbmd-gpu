#include "default_position_controller.h"

#include <thrust/device_ptr.h>

#include "neighbor_list/include/linked_cell/linked_cell_locator.h"
#include "update_position_op.h"
#include <thrust/copy.h>

DefaultPositionController::DefaultPositionController(){};

void DefaultPositionController::Init() 
{
  _dt = DataManager::getInstance().getConfigData()->Get<rbmd::Real>("timestep", "execution");//0.001
}

void DefaultPositionController::Update() {
  bool available_shake = true;
  if (available_shake) {
    bool shake = DataManager::getInstance().getConfigData()->GetJudge<bool>( "fix_shake", "hyper_parameters", "extend");
    if (shake) {
        thrust::copy(_device_data->_d_px.begin(), _device_data->_d_px.end(), _device_data->_d_shake_px.begin());
        thrust::copy(_device_data->_d_py.begin(), _device_data->_d_py.end(), _device_data->_d_shake_py.begin());
        thrust::copy(_device_data->_d_pz.begin(), _device_data->_d_pz.end(), _device_data->_d_shake_pz.begin());

        std::cout<<"---------更新位置内部--------------"<<std::endl;
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
    op::UpdatePositionOp<device::DEVICE_GPU> update_position_op;
    update_position_op(*(_structure_info_data->_num_atoms), _dt, _device_data->_d_box,
                       thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                       thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                       thrust::raw_pointer_cast(_device_data->_d_vz.data()),
                       thrust::raw_pointer_cast(_device_data->_d_px.data()),
                       thrust::raw_pointer_cast(_device_data->_d_py.data()),
                       thrust::raw_pointer_cast(_device_data->_d_pz.data()));
  } else {
    op::UpdatePositionFlagOp<device::DEVICE_GPU> update_position_op;
    update_position_op(*(_structure_info_data->_num_atoms), _dt, _device_data->_d_box,
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
