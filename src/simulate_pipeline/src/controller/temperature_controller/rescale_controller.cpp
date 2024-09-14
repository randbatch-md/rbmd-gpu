#include "rescale_controller.h"

#include <thrust/device_ptr.h>

#include "device_types.h"
#include "neighbor_list/include/linked_cell/linked_cell_locator.h"
#include "rbmd_define.h"
#include "unit_factor.h"
#include "update_temperature_op.h"
#define TEMP_OUTPUT

RescaleController::RescaleController(){};

void RescaleController::Init() {
  _num_atoms = *(_structure_info_data->_num_atoms);
  _temp_sum = 0;

  auto unit = "LJ";                          // �����ļ���ȡ
  UNIT unit_factor = unit_factor_map[unit];  // ����������ض�������

  switch (unit_factor) {
    case UNIT::LJ:
      _mvv2e = UnitFactor<UNIT::LJ>::_mvv2e;
      _kB = UnitFactor<UNIT::LJ>::_kb;
      break;

    case UNIT::REAL:
      _mvv2e = UnitFactor<UNIT::REAL>::_mvv2e;
      _kB = UnitFactor<UNIT::REAL>::_kb;
      break;

    default:
      break;
  }
}

void RescaleController::Update() {
  ComputeTemp();

  UpdataVelocity();
}

void RescaleController::ComputeTemp() {
  rbmd::Real* temp_contrib;

  CHECK_RUNTIME(MALLOC(&temp_contrib, sizeof(rbmd::Real)));
  CHECK_RUNTIME(MEMCPY(temp_contrib, &_temp_sum, sizeof(rbmd::Real), H2D));

  op::ComputeTemperatureOp<device::DEVICE_GPU> compute_temperature_op;
  compute_temperature_op(
      _num_atoms, _mvv2e,
      thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
      thrust::raw_pointer_cast(_device_data->_d_mass.data()),
      thrust::raw_pointer_cast(_device_data->_d_vx.data()),
      thrust::raw_pointer_cast(_device_data->_d_vy.data()),
      thrust::raw_pointer_cast(_device_data->_d_vz.data()), temp_contrib);

  CHECK_RUNTIME(MEMCPY(&_temp_sum, temp_contrib, sizeof(rbmd::Real), D2H));
  // std::cout << "_temp_sum" << _temp_sum << std::endl;

  bool available_shake = false;

  if (available_shake)  // H2O / NACl / EAM ...
  {
    bool shake = true;
    if (shake) {
      _temp = 0.5 * _temp_sum / ((3 * _num_atoms - _num_atoms - 3) * _kB / 2.0);
    } else {
      _temp = 0.5 * _temp_sum / ((3 * _num_atoms - 3) * _kB / 2.0);
    }
  } else  // PEO
  {
    _temp = 0.5 * _temp_sum / ((3 * _num_atoms - 3) * _kB / 2.0);
  }

  std::cout << "_temp=" << _temp << std::endl;
}

void RescaleController::UpdataVelocity() {
  rbmd::Real kbT = 1;  // �����ļ���ȡ
  rbmd::Real coeff_rescale = std::sqrt(kbT / _temp);

  op::UpdataVelocityRescaleOp<device::DEVICE_GPU> updata_velocity_op;
  updata_velocity_op(_num_atoms, coeff_rescale,
                     thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                     thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                     thrust::raw_pointer_cast(_device_data->_d_vz.data()));

#ifdef TEMP_OUTPUT
  // ���������ڴ�
  rbmd::Real* h_vx = new rbmd::Real[_device_data->_d_vx.size()];
  rbmd::Real* h_vy = new rbmd::Real[_device_data->_d_vy.size()];
  rbmd::Real* h_vz = new rbmd::Real[_device_data->_d_vz.size()];

  // ���豸���ݿ���������
  CHECK_RUNTIME(MEMCPY(h_vx,
                       thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                       _device_data->_d_vx.size() * sizeof(rbmd::Real), D2H));
  CHECK_RUNTIME(MEMCPY(h_vy,
                       thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                       _device_data->_d_vy.size() * sizeof(rbmd::Real), D2H));
  CHECK_RUNTIME(MEMCPY(h_vz,
                       thrust::raw_pointer_cast(_device_data->_d_vz.data()),
                       _device_data->_d_vz.size() * sizeof(rbmd::Real), D2H));
  std::ofstream output_file_temp;
  thrust::host_vector<rbmd::Id> atomid2idx =
      LinkedCellLocator::GetInstance().GetLinkedCell()->_atom_id_to_idx;
  output_file_temp.open("output_file_temp_velocity_" +
                        std::to_string(test_current_step) + ".csv");

  // 写入表头
  output_file_temp << "atomid,vx,vy,vz" << std::endl;

  for (size_t i = 0; i < _num_atoms; i++) {
    output_file_temp << i << "," << h_vx[atomid2idx[i]] << "," << h_vy[atomid2idx[i]] << "," << h_vz[atomid2idx[i]] << std::endl;
  }


#endif // V_OUTPUT
}
