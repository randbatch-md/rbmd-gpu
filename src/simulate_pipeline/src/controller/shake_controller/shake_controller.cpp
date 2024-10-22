#include "shake_controller.h"
#include <thrust/device_ptr.h>

#include "data_manager.h"
#include "device_types.h"
#include "neighbor_list/include/linked_cell/linked_cell_locator.h"
#include "unit_factor.h"
#include "shake_controller_op.h"

ShakeController::ShakeController() 
    : _device_data(DataManager::getInstance().getDeviceData())
    , _structure_info_data(DataManager::getInstance().getMDData()->_structure_info_data) {};

void ShakeController::Init()
{
    _num_angle = *(_structure_info_data->_num_angles);

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

void ShakeController::ShakeA() 
{
    op::ShakeAOp<device::DEVICE_GPU> shakeA_op;
    shakeA_op(_num_angle,_dt,_fmt2v,
              _device_data->_d_box,
              (*_structure_info_data->_range)[0][0],
              (*_structure_info_data->_range)[1][0],
              (*_structure_info_data->_range)[2][0],
              (*_structure_info_data->_range)[0][1],
              (*_structure_info_data->_range)[1][1],
              (*_structure_info_data->_range)[2][1],
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
    
    // _device_data->_d_vx = _d_device_data->_shake_vx;
    // _device_data->_d_vy = _d_device_data->_shake_vy;
    // _device_data->_d_vz = _d_device_data->_shake_vz;
    
    // _device_data->_d_px = _d_device_data->_shake_px;
    // _device_data->_d_py = _d_device_data->_shake_py;
    // _device_data->_d_pz = _d_device_data->_shake_pz;
}

void ShakeController::ShakeB() 
{

 // _device_data->_d_vx = _d_device_data->_shake_vx;
 // _device_data->_d_vy = _d_device_data->_shake_vy;
 // _device_data->_d_vz = _d_device_data->_shake_vz;
}