#include "default_position_controller.h"

#include <thrust/device_ptr.h>
#include "unit_factor.h"
#include "data_manager.h"
#include "neighbor_list/include/linked_cell/linked_cell_locator.h"
#include "update_position_op.h"
#include <thrust/copy.h>

DefaultPositionController::DefaultPositionController(){};

void DefaultPositionController::Init() 
{
  auto& num_atoms = *(_structure_info_data->_num_atoms);
  _d_prev_fx.resize(num_atoms, 0.0);
  _d_prev_fy.resize(num_atoms, 0.0);
  _d_prev_fz.resize(num_atoms, 0.0);
  _d_prev_px.resize(num_atoms, 0.0);
  _d_prev_py.resize(num_atoms, 0.0);
  _d_prev_pz.resize(num_atoms, 0.0);
  _dt = DataManager::getInstance().getConfigData()->
    Get<rbmd::Real>("timestep", "execution");//0.001
  _par_a = DataManager::getInstance().getConfigData()->Get<rbmd::Real>(
          "par_a", "execution");
  _par_b = DataManager::getInstance().getConfigData()->Get<rbmd::Real>(
          "par_b", "execution");
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
  _integration_type = DataManager::getInstance().getConfigData()->Get
<std::string>("integration_type", "execution");
  auto d_type_array=DataManager::getInstance().getConfigData()->
  GetArray<rbmd::Real>("d_type", "execution"); //[1.0,1.0,0.1]
  _d1 = d_type_array[0];
  _d2 = d_type_array[1];
  _d3 = d_type_array[2];
  _d4 = d_type_array[3];
  _current_step = 0;
}

void DefaultPositionController::Update() {
  if ("vv" == _integration_type) {
    bool shake = DataManager::getInstance().getConfigData()->GetJudge<bool>
    ( "fix_shake", "hyper_parameters", "extend");
    if (shake) {
        thrust::copy(_device_data->_d_px.begin(), _device_data->_d_px.end(),
          _device_data->_d_shake_px.begin());
        thrust::copy(_device_data->_d_py.begin(), _device_data->_d_py.end(),
          _device_data->_d_shake_py.begin());
        thrust::copy(_device_data->_d_pz.begin(), _device_data->_d_pz.end(),
          _device_data->_d_shake_pz.begin());

      // op::UpdatePositionOp<device::DEVICE_GPU> ()(
      //                   *(_structure_info_data->_num_atoms), _dt,*_box,
      //                    thrust::raw_pointer_cast(_device_data->_d_vx.data()),
      //                    thrust::raw_pointer_cast(_device_data->_d_vy.data()),
      //                    thrust::raw_pointer_cast(_device_data->_d_vz.data()),
      //                    thrust::raw_pointer_cast(_device_data->_d_px.data()),
      //                    thrust::raw_pointer_cast(_device_data->_d_py.data()),
      //                    thrust::raw_pointer_cast(_device_data->_d_pz.data()));
    }
    else {
    op::UpdatePositionFlagOp<device::DEVICE_GPU>()(
                        *(_structure_info_data->_num_atoms), _dt,_fmt2v,*_box,
                       thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
                       thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                       thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                       thrust::raw_pointer_cast(_device_data->_d_vz.data()),
                       thrust::raw_pointer_cast(_device_data->_d_mass.data()),
                       thrust::raw_pointer_cast(_device_data->_d_fx.data()),
                       thrust::raw_pointer_cast(_device_data->_d_fy.data()),
                       thrust::raw_pointer_cast(_device_data->_d_fz.data()),
                       thrust::raw_pointer_cast(_device_data->_d_px.data()),
                       thrust::raw_pointer_cast(_device_data->_d_py.data()),
                       thrust::raw_pointer_cast(_device_data->_d_pz.data()),
                       thrust::raw_pointer_cast(_device_data->_d_flagX.data()),
                       thrust::raw_pointer_cast(_device_data->_d_flagY.data()),
                       thrust::raw_pointer_cast(_device_data->_d_flagZ.data()));
  }
  }
  else if ("test" == _integration_type) {
    op::UpdatePositionFlagOp<device::DEVICE_GPU>()(
                            *(_structure_info_data->_num_atoms), _dt,_fmt2v,*_box,
                           thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
                           thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                           thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                           thrust::raw_pointer_cast(_device_data->_d_vz.data()),
                           thrust::raw_pointer_cast(_device_data->_d_mass.data()),
                           thrust::raw_pointer_cast(_device_data->_d_fx.data()),
                           thrust::raw_pointer_cast(_device_data->_d_fy.data()),
                           thrust::raw_pointer_cast(_device_data->_d_fz.data()),
                           thrust::raw_pointer_cast(_device_data->_d_px.data()),
                           thrust::raw_pointer_cast(_device_data->_d_py.data()),
                           thrust::raw_pointer_cast(_device_data->_d_pz.data()),
                           thrust::raw_pointer_cast(_device_data->_d_flagX.data()),
                           thrust::raw_pointer_cast(_device_data->_d_flagY.data()),
                           thrust::raw_pointer_cast(_device_data->_d_flagZ.data()));

  }
}

void DefaultPositionController::Update1() {
  if ("test" == _integration_type) {
    // std::cout << _d1 << std::endl;
  op::UpdatePositionFlagOp1<device::DEVICE_GPU>()(
                        *(_structure_info_data->_num_atoms),_d1, _dt,*_box,
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

void DefaultPositionController::Update2() {
  if ("test" == _integration_type) {
    op::UpdatePositionFlagOp2<device::DEVICE_GPU>()(
                          *(_structure_info_data->_num_atoms),_d2, _dt,*_box,
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

void DefaultPositionController::Update3() {
  if ("test" == _integration_type) {
    op::UpdatePositionFlagOp3<device::DEVICE_GPU>()(
                          *(_structure_info_data->_num_atoms),_d3,_dt,*_box,
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

void DefaultPositionController::Update4() {
  if ("test" == _integration_type) {
    op::UpdatePositionFlagOp4<device::DEVICE_GPU>()(
                          *(_structure_info_data->_num_atoms),_d4, _dt,*_box,
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

void DefaultPositionController::Updatevl() {
  // _current_step += 1;
  // std::cout << _current_step << std::endl;
    op::UpdatePositionFlagOpvl<device::DEVICE_GPU>()(
         *(_structure_info_data->_num_atoms),_fmt2v,_par_a,_par_b, _dt,test_current_step,*_box,
         thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
         thrust::raw_pointer_cast(_device_data->_d_fx.data()),
         thrust::raw_pointer_cast(_device_data->_d_fy.data()),
         thrust::raw_pointer_cast(_device_data->_d_fz.data()),
         thrust::raw_pointer_cast(_device_data->_d_mass.data()),
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

void DefaultPositionController::Updatebm() {
  // _current_step += 1;
  // std::cout << _current_step << std::endl;
  op::UpdatePositionFlagOpbm<device::DEVICE_GPU>()(
       *(_structure_info_data->_num_atoms),_fmt2v, _par_a,_par_b, _dt, test_current_step,*_box,
       thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
       thrust::raw_pointer_cast(_device_data->_d_fx.data()),
       thrust::raw_pointer_cast(_device_data->_d_fy.data()),
       thrust::raw_pointer_cast(_device_data->_d_fz.data()),
       thrust::raw_pointer_cast(_device_data->_d_mass.data()),
       thrust::raw_pointer_cast(_d_prev_fx.data()),
       thrust::raw_pointer_cast(_d_prev_fy.data()),
       thrust::raw_pointer_cast(_d_prev_fz.data()),
       thrust::raw_pointer_cast(_device_data->_d_vx.data()),
       thrust::raw_pointer_cast(_device_data->_d_vy.data()),
       thrust::raw_pointer_cast(_device_data->_d_vz.data()),
       thrust::raw_pointer_cast(_device_data->_d_px.data()),
       thrust::raw_pointer_cast(_device_data->_d_py.data()),
       thrust::raw_pointer_cast(_device_data->_d_pz.data()),
       thrust::raw_pointer_cast(_device_data->_d_flagX.data()),
       thrust::raw_pointer_cast(_device_data->_d_flagY.data()),
       thrust::raw_pointer_cast(_device_data->_d_flagZ.data()));
  // 同步并打印设备数据
  // CHECK_RUNTIME(DEVICESYNC());
  // thrust::host_vector<rbmd::Real> h_fx_prev = _d_prev_fx;
  // printf("CPU: h_fx_prev[1] = %f\n", h_fx_prev[1]);

}


void DefaultPositionController::SetCenterTargetPositions() {
  std::string init_type = "inbuild";
  if (init_type == _init_type) {
  }
}
