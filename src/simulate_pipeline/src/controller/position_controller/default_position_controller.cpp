#include "default_position_controller.h"

#include <thrust/copy.h>
#include <thrust/device_ptr.h>

#include "group_controller_op.h"
#include "neighbor_list/include/linked_cell/linked_cell_locator.h"
#include "simulate.h"
#include "unit_factor.h"
#include "update_position_op.h"

DefaultPositionController::DefaultPositionController(){};

void DefaultPositionController::Init() 
{
  auto& num_atoms = *(_structure_info_data->_num_atoms);
  _dt = DataManager::getInstance().getConfigData()->
    Get<rbmd::Real>("timestep", "execution");//0.001

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

  const auto& config = DataManager::getInstance().getConfigData();
  auto integration_type = config->Get<std::string>("integration_type", "execution");
  if("bm" ==integration_type){
    _par_a = DataManager::getInstance().getConfigData()->Get<rbmd::Real>(
        "par_a", "execution");
    _par_b = DataManager::getInstance().getConfigData()->Get<rbmd::Real>(
            "par_b", "execution");
    _d_prev_fx.resize(num_atoms, 0.0);
    _d_prev_fy.resize(num_atoms, 0.0);
    _d_prev_fz.resize(num_atoms, 0.0);
    _d_prev_px.resize(num_atoms, 0.0);
    _d_prev_py.resize(num_atoms, 0.0);
    _d_prev_pz.resize(num_atoms, 0.0);
  }
}

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
                        *(_structure_info_data->_num_atoms), _dt,*_box,
                       thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                       thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                       thrust::raw_pointer_cast(_device_data->_d_vz.data()),
                       thrust::raw_pointer_cast(_device_data->_d_px.data()),
                       thrust::raw_pointer_cast(_device_data->_d_py.data()),
                       thrust::raw_pointer_cast(_device_data->_d_pz.data()),
                       thrust::raw_pointer_cast(_device_data->_d_flagX.data()),
                       thrust::raw_pointer_cast(_device_data->_d_flagY.data()),
                       thrust::raw_pointer_cast(_device_data->_d_flagZ.data()));

      //Unwarp Position
      op::UnwarpPositionOp<device::DEVICE_GPU>()(*(_structure_info_data->_num_atoms),*_box,
        thrust::raw_pointer_cast(_device_data->_d_px.data()),
        thrust::raw_pointer_cast(_device_data->_d_py.data()),
        thrust::raw_pointer_cast(_device_data->_d_pz.data()),
        thrust::raw_pointer_cast(_device_data->_d_flagX.data()),
        thrust::raw_pointer_cast(_device_data->_d_flagY.data()),
        thrust::raw_pointer_cast(_device_data->_d_flagZ.data()),
        thrust::raw_pointer_cast(_device_data->_d_unwarp_px.data()),
        thrust::raw_pointer_cast(_device_data->_d_unwarp_py.data()),
        thrust::raw_pointer_cast(_device_data->_d_unwarp_pz.data()));
  }
}

void DefaultPositionController::Updatebm() {
  // _current_step += 1;
  //std::cout << "test_current_step--p: "  << test_current_step <<  std::endl;
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
}

void DefaultPositionController::SetCenterTargetPositions() {
  std::string init_type = "inbuild";
  if (init_type == _init_type) {
  }
}

void  DefaultPositionController::PBC() {
  op::PBCOp<device::DEVICE_GPU>()(
    *(_structure_info_data->_num_atoms),*_box,
    thrust::raw_pointer_cast(_device_data->_d_px.data()),
    thrust::raw_pointer_cast(_device_data->_d_py.data()),
    thrust::raw_pointer_cast(_device_data->_d_pz.data()),
    thrust::raw_pointer_cast(_device_data->_d_flagX.data()),
    thrust::raw_pointer_cast(_device_data->_d_flagY.data()),
    thrust::raw_pointer_cast(_device_data->_d_flagZ.data()));
}