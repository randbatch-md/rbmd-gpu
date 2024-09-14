#include "default_velocity_controller.h"

#include <thrust/device_ptr.h>

#include "data_manager.h"
#include "device_types.h"
#include "neighbor_list/include/linked_cell/linked_cell_locator.h"
#include "unit_factor.h"
#include "update_velocity_op.h"
#define V_OUTPUT

DefaultVelocityController::DefaultVelocityController(){};

void DefaultVelocityController::Init() {
  _num_atoms = *(_structure_info_data->_num_atoms);

  _dt = 0.001;
  auto unit = "LJ";
  UNIT unit_factor = unit_factor_map[unit];

  switch (unit_factor) {
    case UNIT::METAL:
      _fmt2v = UnitFactor<UNIT::METAL>::_kb;
      break;
    case UNIT::LJ:
      _fmt2v = UnitFactor<UNIT::LJ>::_kb;
      break;
    case UNIT::REAL:
      _fmt2v = UnitFactor<UNIT::REAL>::_kb;
      break;
    default:
      break;
  }
}

void DefaultVelocityController::Update() {
  bool shake = false;  
  if (shake) {
    // __device_data->_shake_vx = _device_data->_d_vx;
    // __device_data->_shake_vy = _device_data->_d_vy;
    // __device_data->_shake_vz = _device_data->_d_vz;
  }

  op::UpdateVelocityOp<device::DEVICE_GPU> update_velocity_op;
  update_velocity_op(_num_atoms, _dt, _fmt2v,
                     thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
                     thrust::raw_pointer_cast(_device_data->_d_mass.data()),
                     thrust::raw_pointer_cast(_device_data->_d_fx.data()),
                     thrust::raw_pointer_cast(_device_data->_d_fy.data()),
                     thrust::raw_pointer_cast(_device_data->_d_fz.data()),
                     thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                     thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                     thrust::raw_pointer_cast(_device_data->_d_vz.data()));
  
#ifdef V_OUTPUT
  rbmd::Real* h_vx = new rbmd::Real[_device_data->_d_vx.size()];
  rbmd::Real* h_vy = new rbmd::Real[_device_data->_d_vy.size()];
  rbmd::Real* h_vz = new rbmd::Real[_device_data->_d_vz.size()];

  rbmd::Real* h_fx = new rbmd::Real[_device_data->_d_fx.size()];
  rbmd::Real* h_fy = new rbmd::Real[_device_data->_d_fy.size()];
  rbmd::Real* h_fz = new rbmd::Real[_device_data->_d_fz.size()];

  CHECK_RUNTIME(MEMCPY(h_vx, thrust::raw_pointer_cast(_device_data->_d_vx.data()), _device_data->_d_vx.size() * sizeof(rbmd::Real), D2H));
  CHECK_RUNTIME(MEMCPY(h_vy, thrust::raw_pointer_cast(_device_data->_d_vy.data()), _device_data->_d_vy.size() * sizeof(rbmd::Real), D2H));
  CHECK_RUNTIME(MEMCPY(h_vz, thrust::raw_pointer_cast(_device_data->_d_vz.data()), _device_data->_d_vz.size() * sizeof(rbmd::Real), D2H));

  CHECK_RUNTIME(MEMCPY(h_fx, thrust::raw_pointer_cast(_device_data->_d_fx.data()), _device_data->_d_fx.size() * sizeof(rbmd::Real), D2H));
  CHECK_RUNTIME(MEMCPY(h_fy, thrust::raw_pointer_cast(_device_data->_d_fy.data()), _device_data->_d_fy.size() * sizeof(rbmd::Real), D2H));
  CHECK_RUNTIME(MEMCPY(h_fz, thrust::raw_pointer_cast(_device_data->_d_fz.data()), _device_data->_d_fz.size() * sizeof(rbmd::Real), D2H));

  std::ofstream output_file_v1_velocity;
  thrust::host_vector<rbmd::Id> atomid2idx=  LinkedCellLocator::GetInstance().GetLinkedCell()->_atom_id_to_idx;
  output_file_v1_velocity.open("output_file_v1_velocity_" + std::to_string(test_current_step) + ".csv");

  // 写入表头
  output_file_v1_velocity << "atomid,vx,vy,vz" << std::endl;

  for (size_t i = 0; i < _num_atoms; i++)
  {
    output_file_v1_velocity << i << "," << h_vx[atomid2idx[i]] << "," << h_vy[atomid2idx[i]] << "," << h_vz[atomid2idx[i]] << std::endl;
  }

  std::ofstream output_file_v1_force;
  output_file_v1_force.open("output_file_v1_force_" + std::to_string(test_current_step) + ".csv");

  // 写入表头
  output_file_v1_force << "atomid,fx,fy,fz" << std::endl;

  for (size_t i = 0; i < _num_atoms; i++)
  {
    output_file_v1_force << i << "," << h_fx[atomid2idx[i]] << "," << h_fy[atomid2idx[i]] << "," << h_fz[atomid2idx[i]] << std::endl;
  }

#endif // V_OUTPUT


}

void DefaultVelocityController::Update2() {
    bool shake = false;
    if (shake) {
        // __device_data->_shake_vx = _device_data->_d_vx;
        // __device_data->_shake_vy = _device_data->_d_vy;
        // __device_data->_shake_vz = _device_data->_d_vz;
    }

    op::UpdateVelocityOp<device::DEVICE_GPU> update_velocity_op;
    update_velocity_op(_num_atoms, _dt, _fmt2v,
        thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
        thrust::raw_pointer_cast(_device_data->_d_mass.data()),
        thrust::raw_pointer_cast(_device_data->_d_fx.data()),
        thrust::raw_pointer_cast(_device_data->_d_fy.data()),
        thrust::raw_pointer_cast(_device_data->_d_fz.data()),
        thrust::raw_pointer_cast(_device_data->_d_vx.data()),
        thrust::raw_pointer_cast(_device_data->_d_vy.data()),
        thrust::raw_pointer_cast(_device_data->_d_vz.data()));

#ifdef V_OUTPUT
    rbmd::Real* h_vx = new rbmd::Real[_device_data->_d_vx.size()];
    rbmd::Real* h_vy = new rbmd::Real[_device_data->_d_vy.size()];
    rbmd::Real* h_vz = new rbmd::Real[_device_data->_d_vz.size()];

    rbmd::Real* h_fx = new rbmd::Real[_device_data->_d_fx.size()];
    rbmd::Real* h_fy = new rbmd::Real[_device_data->_d_fy.size()];
    rbmd::Real* h_fz = new rbmd::Real[_device_data->_d_fz.size()];

    CHECK_RUNTIME(MEMCPY(h_vx, thrust::raw_pointer_cast(_device_data->_d_vx.data()), _device_data->_d_vx.size() * sizeof(rbmd::Real), D2H));
    CHECK_RUNTIME(MEMCPY(h_vy, thrust::raw_pointer_cast(_device_data->_d_vy.data()), _device_data->_d_vy.size() * sizeof(rbmd::Real), D2H));
    CHECK_RUNTIME(MEMCPY(h_vz, thrust::raw_pointer_cast(_device_data->_d_vz.data()), _device_data->_d_vz.size() * sizeof(rbmd::Real), D2H));

    CHECK_RUNTIME(MEMCPY(h_fx, thrust::raw_pointer_cast(_device_data->_d_fx.data()), _device_data->_d_fx.size() * sizeof(rbmd::Real), D2H));
    CHECK_RUNTIME(MEMCPY(h_fy, thrust::raw_pointer_cast(_device_data->_d_fy.data()), _device_data->_d_fy.size() * sizeof(rbmd::Real), D2H));
    CHECK_RUNTIME(MEMCPY(h_fz, thrust::raw_pointer_cast(_device_data->_d_fz.data()), _device_data->_d_fz.size() * sizeof(rbmd::Real), D2H));

  std::ofstream output_file_v2_velocity;
  thrust::host_vector<rbmd::Id> atomid2idx=  LinkedCellLocator::GetInstance().GetLinkedCell()->_atom_id_to_idx;
  output_file_v2_velocity.open("output_file_v1_velocity_" + std::to_string(test_current_step) + ".csv");

  // 写入表头
  output_file_v2_velocity << "atomid,vx,vy,vz" << std::endl;

  for (size_t i = 0; i < _num_atoms; i++)
  {
    output_file_v2_velocity << i << "," << h_vx[atomid2idx[i]] << "," << h_vy[atomid2idx[i]] << "," << h_vz[atomid2idx[i]] << std::endl;
  }

  std::ofstream output_file_v2_force;
  output_file_v2_force.open("output_file_v1_force_" + std::to_string(test_current_step) + ".csv");

  // 写入表头
  output_file_v2_force << "atomid,fx,fy,fz" << std::endl;

  for (size_t i = 0; i < _num_atoms; i++)
  {
    output_file_v2_force << i << "," << h_fx[atomid2idx[i]] << "," << h_fy[atomid2idx[i]] << "," << h_fz[atomid2idx[i]] << std::endl;
  }

#endif // V_OUTPUT


}