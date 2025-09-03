#include "default_position_controller.h"

#include <thrust/device_ptr.h>

#include "data_manager.h"
#include "neighbor_list/include/linked_cell/linked_cell_locator.h"
#include "unit_factor.h"
#include "update_position_op.h"
#include <thrust/copy.h>

DefaultPositionController::DefaultPositionController(){};

void DefaultPositionController::Init() 
{
  _dt = DataManager::getInstance().getConfigData()->
    Get<rbmd::Real>("timestep", "execution");//0.001

  _integration_type = DataManager::getInstance().getConfigData()->Get
<std::string>("integration_type", "execution");
  auto d_type_array=DataManager::getInstance().getConfigData()->
  GetArray<rbmd::Real>("d_type", "execution"); //[1.0,1.0,0.1]
  _d1 = d_type_array[0];
  _d2 = d_type_array[1];
  _d3 = d_type_array[2];
  _d4 = d_type_array[3];

  // ... 您已有的代码 ...
  auto num_atoms = *(_structure_info_data->_num_atoms);

  // VVVV 在这里添加 VVVV
  // 为这个控制器私有的历史向量分配内存
  _d_prev_fx.resize(num_atoms,0.0);
  _d_prev_fy.resize(num_atoms,0.0);
  _d_prev_fz.resize(num_atoms,0.0);
  // ^^^^ 添加结束 ^^^^
}

//leapfrog
void DefaultPositionController::Update() {
    bool shake = DataManager::getInstance().getConfigData()->GetJudge<bool>
( "fix_shake", "hyper_parameters", "extend");
    if (shake) {
      thrust::copy(_device_data->_d_px.begin(), _device_data->_d_px.end(),
        _device_data->_d_shake_px.begin());
      thrust::copy(_device_data->_d_py.begin(), _device_data->_d_py.end(),
        _device_data->_d_shake_py.begin());
      thrust::copy(_device_data->_d_pz.begin(), _device_data->_d_pz.end(),
        _device_data->_d_shake_pz.begin());
      op::UpdatePositionOp<device::DEVICE_GPU> ()(
                        *(_structure_info_data->_num_atoms), _dt,*_box,
                         thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                         thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                         thrust::raw_pointer_cast(_device_data->_d_vz.data()),
                         thrust::raw_pointer_cast(_device_data->_d_px.data()),
                         thrust::raw_pointer_cast(_device_data->_d_py.data()),
                         thrust::raw_pointer_cast(_device_data->_d_pz.data()));
    }
    else {
      op::UpdatePositionFlagOp<device::DEVICE_GPU>()(
                          *(_structure_info_data->_num_atoms), _dt, _fmt2v,
                         thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
                         thrust::raw_pointer_cast(_device_data->_d_mass.data()),
                          *_box,
                         thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                         thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                         thrust::raw_pointer_cast(_device_data->_d_vz.data()),
                         thrust::raw_pointer_cast(_device_data->_d_px.data()),
                         thrust::raw_pointer_cast(_device_data->_d_py.data()),
                         thrust::raw_pointer_cast(_device_data->_d_pz.data()),
                         thrust::raw_pointer_cast(_device_data->_d_flagX.data()),
                         thrust::raw_pointer_cast(_device_data->_d_flagY.data()),
                         thrust::raw_pointer_cast(_device_data->_d_flagZ.data()),
                         thrust::raw_pointer_cast(_device_data->_d_fx.data()),
                         thrust::raw_pointer_cast(_device_data->_d_fy.data()),
                         thrust::raw_pointer_cast(_device_data->_d_fz.data()));
    }
  }

//vv
void DefaultPositionController::Update_vv() {
  op::UpdatePositionFlagOpvv<device::DEVICE_GPU>()(
                      *(_structure_info_data->_num_atoms), _dt, _fmt2v,
                     thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
                     thrust::raw_pointer_cast(_device_data->_d_mass.data()),
                      *_box,
                     thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                     thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                     thrust::raw_pointer_cast(_device_data->_d_vz.data()),
                     thrust::raw_pointer_cast(_device_data->_d_px.data()),
                     thrust::raw_pointer_cast(_device_data->_d_py.data()),
                     thrust::raw_pointer_cast(_device_data->_d_pz.data()),
                     thrust::raw_pointer_cast(_device_data->_d_flagX.data()),
                     thrust::raw_pointer_cast(_device_data->_d_flagY.data()),
                     thrust::raw_pointer_cast(_device_data->_d_flagZ.data()),
                     thrust::raw_pointer_cast(_device_data->_d_fx.data()),
                     thrust::raw_pointer_cast(_device_data->_d_fy.data()),
                     thrust::raw_pointer_cast(_device_data->_d_fz.data()));
}

// PRK
void DefaultPositionController::Update1() {
    op::UpdatePositionFlagOp1<device::DEVICE_GPU>()(
                        *(_structure_info_data->_num_atoms), _d1, _dt, _fmt2v,
                       thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
                       thrust::raw_pointer_cast(_device_data->_d_mass.data()),
                        *_box,
                       thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                       thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                       thrust::raw_pointer_cast(_device_data->_d_vz.data()),
                       thrust::raw_pointer_cast(_device_data->_d_px.data()),
                       thrust::raw_pointer_cast(_device_data->_d_py.data()),
                       thrust::raw_pointer_cast(_device_data->_d_pz.data()),
                       thrust::raw_pointer_cast(_device_data->_d_flagX.data()),
                       thrust::raw_pointer_cast(_device_data->_d_flagY.data()),
                       thrust::raw_pointer_cast(_device_data->_d_flagZ.data()),
                       thrust::raw_pointer_cast(_device_data->_d_fx.data()),
                       thrust::raw_pointer_cast(_device_data->_d_fy.data()),
                       thrust::raw_pointer_cast(_device_data->_d_fz.data()));
}
void DefaultPositionController::Update2() {
    op::UpdatePositionFlagOp2<device::DEVICE_GPU>()(
                        *(_structure_info_data->_num_atoms), _d2, _dt, _fmt2v,
                       thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
                       thrust::raw_pointer_cast(_device_data->_d_mass.data()),
                        *_box,
                       thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                       thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                       thrust::raw_pointer_cast(_device_data->_d_vz.data()),
                       thrust::raw_pointer_cast(_device_data->_d_px.data()),
                       thrust::raw_pointer_cast(_device_data->_d_py.data()),
                       thrust::raw_pointer_cast(_device_data->_d_pz.data()),
                       thrust::raw_pointer_cast(_device_data->_d_flagX.data()),
                       thrust::raw_pointer_cast(_device_data->_d_flagY.data()),
                       thrust::raw_pointer_cast(_device_data->_d_flagZ.data()),
                       thrust::raw_pointer_cast(_device_data->_d_fx.data()),
                       thrust::raw_pointer_cast(_device_data->_d_fy.data()),
                       thrust::raw_pointer_cast(_device_data->_d_fz.data()));
}
void DefaultPositionController::Update3() {
    op::UpdatePositionFlagOp3<device::DEVICE_GPU>()(
                        *(_structure_info_data->_num_atoms), _d3, _dt, _fmt2v,
                       thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
                       thrust::raw_pointer_cast(_device_data->_d_mass.data()),
                        *_box,
                       thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                       thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                       thrust::raw_pointer_cast(_device_data->_d_vz.data()),
                       thrust::raw_pointer_cast(_device_data->_d_px.data()),
                       thrust::raw_pointer_cast(_device_data->_d_py.data()),
                       thrust::raw_pointer_cast(_device_data->_d_pz.data()),
                       thrust::raw_pointer_cast(_device_data->_d_flagX.data()),
                       thrust::raw_pointer_cast(_device_data->_d_flagY.data()),
                       thrust::raw_pointer_cast(_device_data->_d_flagZ.data()),
                       thrust::raw_pointer_cast(_device_data->_d_fx.data()),
                       thrust::raw_pointer_cast(_device_data->_d_fy.data()),
                       thrust::raw_pointer_cast(_device_data->_d_fz.data()));
}
void DefaultPositionController::Update4() {
    op::UpdatePositionFlagOp4<device::DEVICE_GPU>()(
                        *(_structure_info_data->_num_atoms), _d4, _dt, _fmt2v,
                       thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
                       thrust::raw_pointer_cast(_device_data->_d_mass.data()),
                        *_box,
                       thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                       thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                       thrust::raw_pointer_cast(_device_data->_d_vz.data()),
                       thrust::raw_pointer_cast(_device_data->_d_px.data()),
                       thrust::raw_pointer_cast(_device_data->_d_py.data()),
                       thrust::raw_pointer_cast(_device_data->_d_pz.data()),
                       thrust::raw_pointer_cast(_device_data->_d_flagX.data()),
                       thrust::raw_pointer_cast(_device_data->_d_flagY.data()),
                       thrust::raw_pointer_cast(_device_data->_d_flagZ.data()),
                       thrust::raw_pointer_cast(_device_data->_d_fx.data()),
                       thrust::raw_pointer_cast(_device_data->_d_fy.data()),
                       thrust::raw_pointer_cast(_device_data->_d_fz.data()));
}

void DefaultPositionController::Update_Beeman() {
  op::UpdatePositionFlagOpBeeman<device::DEVICE_GPU>()(
                      *(_structure_info_data->_num_atoms),_dt, test_current_step,_fmt2v,
                     thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
                     thrust::raw_pointer_cast(_device_data->_d_mass.data()),
                      *_box,
                     thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                     thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                     thrust::raw_pointer_cast(_device_data->_d_vz.data()),
                     thrust::raw_pointer_cast(_device_data->_d_fx.data()),
                     thrust::raw_pointer_cast(_device_data->_d_fy.data()),
                     thrust::raw_pointer_cast(_device_data->_d_fz.data()),
                     thrust::raw_pointer_cast(_d_prev_fx.data()),
                     thrust::raw_pointer_cast(_d_prev_fy.data()),
                     thrust::raw_pointer_cast(_d_prev_fz.data()),
                     thrust::raw_pointer_cast(_device_data->_d_px.data()),
                     thrust::raw_pointer_cast(_device_data->_d_py.data()),
                     thrust::raw_pointer_cast(_device_data->_d_pz.data()),
                     thrust::raw_pointer_cast(_device_data->_d_flagX.data()),
                     thrust::raw_pointer_cast(_device_data->_d_flagY.data()),
                     thrust::raw_pointer_cast(_device_data->_d_flagZ.data()));
}

void DefaultPositionController::SetCenterTargetPositions() {
  std::string init_type = "inbuild";
  if (init_type == _init_type) {
  }
}
