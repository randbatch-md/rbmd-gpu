#include "default_velocity_controller.h"

#include <thrust/device_ptr.h>

#include "common/unit_factor.h"
#include "velocity_controller_op/update_velocity_op.h"

DefaultVelocityController::DefaultVelocityController(){};

void DefaultVelocityController::Init() {
  _num_atoms = _structure_info_data->_num_atoms;

  // �����ļ��ж�ȡ
  _dt = 0.001;
  auto unit = "LJ";
  UNIT unit_factor = unit_factor_map[unit];

  switch (unit_factor_map["lj"]) {
    case UNIT::LJ:
      _fmt2v = UnitFactor<UNIT::LJ>::_kb;
    case UNIT::REAL:
      _fmt2v = UnitFactor<UNIT::REAL>::_kb;
    default:
      break;
  }
}

void DefaultVelocityController::Update() {
  bool shake = false;  // �����ļ��ж�ȡ
  if (shake) {
    // ��ǰLJ����Ҫ���shake���������
    // __device_data->_shake_vx = _device_data->_d_vx;
    // __device_data->_shake_vy = _device_data->_d_vy;
    // __device_data->_shake_vz = _device_data->_d_vz;
  }

  op::UpdateVelocityOp<device::DEVICE_GPU> update_velocity_op;
  update_velocity_op(_num_atoms, _dt, _fmt2v,
                     thrust::raw_pointer_cast(_device_data->_d_mass.data()),
                     thrust::raw_pointer_cast(_device_data->_d_fx.data()),
                     thrust::raw_pointer_cast(_device_data->_d_fy.data()),
                     thrust::raw_pointer_cast(_device_data->_d_fz.data()),
                     thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                     thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                     thrust::raw_pointer_cast(_device_data->_d_vz.data()));
}

