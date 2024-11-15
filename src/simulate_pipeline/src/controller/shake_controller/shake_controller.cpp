#include "shake_controller.h"
#include <thrust/device_ptr.h>

#include "data_manager.h"
#include "device_types.h"
#include "neighbor_list/include/linked_cell/linked_cell_locator.h"
#include "unit_factor.h"
#include "shake_controller_op.h"
#include <thrust/copy.h>

ShakeController::ShakeController() 
    : _device_data(DataManager::getInstance().getDeviceData())
    , _structure_info_data(DataManager::getInstance().getMDData()->_structure_info_data) {};

void ShakeController::Init()
{
    _dt = DataManager::getInstance().getConfigData()->Get<rbmd::Real>
    ( "timestep", "execution"); // TODO: Json file
    auto unit = DataManager::getInstance().getConfigData()->Get<std::string>
    ( "unit", "init_configuration", "read_data");
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

void ShakeController::ShakeA() 
{

    auto atom_id_to_idx = LinkedCellLocator::GetInstance().GetLinkedCell()->_atom_id_to_idx;
    op::ShakeAOp<device::DEVICE_GPU> shakeA_op;
    shakeA_op(*(_structure_info_data->_num_angles),_dt,_fmt2v,
              _device_data->_d_box,
              thrust::raw_pointer_cast(atom_id_to_idx.data()),
              thrust::raw_pointer_cast(_device_data->_d_mass.data()),
              thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
              thrust::raw_pointer_cast(_device_data->_d_angle_id_vec.data()),//TODO: qw:Real3 is ok?
              thrust::raw_pointer_cast(_device_data->_d_shake_px.data()),
              thrust::raw_pointer_cast(_device_data->_d_shake_py.data()),
              thrust::raw_pointer_cast(_device_data->_d_shake_pz.data()),
              thrust::raw_pointer_cast(_device_data->_d_shake_vx.data()),
              thrust::raw_pointer_cast(_device_data->_d_shake_vy.data()),
              thrust::raw_pointer_cast(_device_data->_d_shake_vz.data()),
              thrust::raw_pointer_cast(_device_data->_d_fx.data()),
              thrust::raw_pointer_cast(_device_data->_d_fy.data()),
              thrust::raw_pointer_cast(_device_data->_d_fz.data()),
              thrust::raw_pointer_cast(_device_data->_d_flagX.data()),
              thrust::raw_pointer_cast(_device_data->_d_flagY.data()),
              thrust::raw_pointer_cast(_device_data->_d_flagZ.data())); // TODO:FLAG & locator


    thrust::copy(_device_data->_d_shake_vx.begin(), _device_data->_d_shake_vx.end(),
      _device_data->_d_vx.begin());
    thrust::copy(_device_data->_d_shake_vy.begin(), _device_data->_d_shake_vy.end(),
      _device_data->_d_vy.begin());
    thrust::copy(_device_data->_d_shake_vz.begin(), _device_data->_d_shake_vz.end(),
      _device_data->_d_vz.begin());

    thrust::copy(_device_data->_d_shake_px.begin(), _device_data->_d_shake_px.end(),
      _device_data->_d_px.begin());
    thrust::copy(_device_data->_d_shake_py.begin(), _device_data->_d_shake_py.end(),
      _device_data->_d_py.begin());
    thrust::copy(_device_data->_d_shake_pz.begin(), _device_data->_d_shake_pz.end(),
      _device_data->_d_pz.begin());
}

void ShakeController::ShakeB() 
{
    auto atom_id_to_idx = LinkedCellLocator::GetInstance().GetLinkedCell()->_atom_id_to_idx;
    op::ShakeBOp<device::DEVICE_GPU> shakeB_op;
    shakeB_op(*(_structure_info_data->_num_angles),_dt,_fmt2v,
              _device_data->_d_box,
              thrust::raw_pointer_cast(atom_id_to_idx.data()),
              thrust::raw_pointer_cast(_device_data->_d_mass.data()),
              thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
              thrust::raw_pointer_cast(_device_data->_d_angle_id_vec.data()),//TODO: qw:Real3 is ok?
              thrust::raw_pointer_cast(_device_data->_d_px.data()),
              thrust::raw_pointer_cast(_device_data->_d_py.data()),
              thrust::raw_pointer_cast(_device_data->_d_pz.data()),
              thrust::raw_pointer_cast(_device_data->_d_shake_vx.data()),
              thrust::raw_pointer_cast(_device_data->_d_shake_vy.data()),
              thrust::raw_pointer_cast(_device_data->_d_shake_vz.data()),
              thrust::raw_pointer_cast(_device_data->_d_fx.data()),
              thrust::raw_pointer_cast(_device_data->_d_fy.data()),
              thrust::raw_pointer_cast(_device_data->_d_fz.data())); // TODO:FLAG & locator

    thrust::copy(_device_data->_d_shake_vx.begin(), _device_data->_d_shake_vx.end(),
      _device_data->_d_vx.begin());
    thrust::copy(_device_data->_d_shake_vy.begin(), _device_data->_d_shake_vy.end(),
      _device_data->_d_vy.begin());
    thrust::copy(_device_data->_d_shake_vz.begin(), _device_data->_d_shake_vz.end(),
      _device_data->_d_vz.begin());
}