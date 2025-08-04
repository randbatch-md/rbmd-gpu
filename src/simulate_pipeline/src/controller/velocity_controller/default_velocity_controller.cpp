#include "default_velocity_controller.h"

#include <thrust/device_ptr.h>

#include "data_manager.h"
#include "device_types.h"
#include "neighbor_list/include/linked_cell/linked_cell_locator.h"
#include "unit_factor.h"
#include "update_velocity_op.h"
#include <thrust/copy.h>

DefaultVelocityController::DefaultVelocityController(){};

void DefaultVelocityController::Init() {

  _dt = DataManager::getInstance().getConfigData()->Get<rbmd::Real>(
          "timestep", "execution");//0.001
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
  auto c_type_array=DataManager::getInstance().getConfigData()->
  GetArray<rbmd::Real>("c_type", "execution"); //[1.0,1.0,0.1]
  _c1 = c_type_array[0];
  _c2 = c_type_array[1];
  _c3 = c_type_array[2];
  _c4 = c_type_array[3];
}

//leapfrog
void DefaultVelocityController::Update() {
    bool shake = DataManager::getInstance().getConfigData()->GetJudge<bool>
    ( "fix_shake", "hyper_parameters", "extend");
    if (shake) {
      thrust::copy(_device_data->_d_vx.begin(),_device_data->_d_vx.end(),
        _device_data->_d_shake_vx.begin());
      thrust::copy(_device_data->_d_vy.begin(),_device_data->_d_vy.end(),
        _device_data->_d_shake_vy.begin());
      thrust::copy(_device_data->_d_vz.begin(),_device_data->_d_vz.end(),
        _device_data->_d_shake_vz.begin());
    }
    else {
      op::UpdateVelocityOp<device::DEVICE_GPU>()
      (
          *(_structure_info_data->_num_atoms), _dt, _fmt2v,
          thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
          thrust::raw_pointer_cast(_device_data->_d_mass.data()),
          thrust::raw_pointer_cast(_device_data->_d_fx.data()),
          thrust::raw_pointer_cast(_device_data->_d_fy.data()),
          thrust::raw_pointer_cast(_device_data->_d_fz.data()),
          thrust::raw_pointer_cast(_device_data->_d_vx.data()),
          thrust::raw_pointer_cast(_device_data->_d_vy.data()),
          thrust::raw_pointer_cast(_device_data->_d_vz.data()));
    }
  }

//vv
void DefaultVelocityController::Update_vv() {
  // // 获取力的设备指针
  // auto* d_fx = thrust::raw_pointer_cast(_device_data->_d_fx.data());
  // auto* d_fy = thrust::raw_pointer_cast(_device_data->_d_fy.data());
  // auto* d_fz = thrust::raw_pointer_cast(_device_data->_d_fz.data());
  //
  // // 复制前N个力的值到主机内存（例如前5个原子）
  // const int N = 5;
  // std::vector<rbmd::Real> h_fx(N), h_fy(N), h_fz(N);
  //
  // cudaMemcpy(h_fx.data(), d_fx, sizeof(rbmd::Real) * N, cudaMemcpyDeviceToHost);
  // cudaMemcpy(h_fy.data(), d_fy, sizeof(rbmd::Real) * N, cudaMemcpyDeviceToHost);
  // cudaMemcpy(h_fz.data(), d_fz, sizeof(rbmd::Real) * N, cudaMemcpyDeviceToHost);
  //
  // // 打印结果
  // std::cout << "First " << N << " force values:" << std::endl;
  // for (int i = 0; i < N; ++i) {
  //   std::cout << "Atom " << i << ": fx=" << h_fx[i]
  //             << ", fy=" << h_fy[i]
  //             << ", fz=" << h_fz[i] << std::endl;
  // }
  op::UpdateVelocityOpvv<device::DEVICE_GPU>()
  (
      *(_structure_info_data->_num_atoms), _dt, _fmt2v,
      thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
      thrust::raw_pointer_cast(_device_data->_d_mass.data()),
      thrust::raw_pointer_cast(_device_data->_d_fx.data()),
      thrust::raw_pointer_cast(_device_data->_d_fy.data()),
      thrust::raw_pointer_cast(_device_data->_d_fz.data()),
      thrust::raw_pointer_cast(_device_data->_d_vx.data()),
      thrust::raw_pointer_cast(_device_data->_d_vy.data()),
      thrust::raw_pointer_cast(_device_data->_d_vz.data()));
}

//PRK
void DefaultVelocityController::Update1() {
    op::UpdateVelocityOp1<device::DEVICE_GPU>()
    (
        *(_structure_info_data->_num_atoms), _c1, _dt, _fmt2v,
        thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
        thrust::raw_pointer_cast(_device_data->_d_mass.data()),
        thrust::raw_pointer_cast(_device_data->_d_fx.data()),
        thrust::raw_pointer_cast(_device_data->_d_fy.data()),
        thrust::raw_pointer_cast(_device_data->_d_fz.data()),
        thrust::raw_pointer_cast(_device_data->_d_vx.data()),
        thrust::raw_pointer_cast(_device_data->_d_vy.data()),
        thrust::raw_pointer_cast(_device_data->_d_vz.data()));
}

void DefaultVelocityController::Update2() {
    op::UpdateVelocityOp2<device::DEVICE_GPU>()
    (
        *(_structure_info_data->_num_atoms), _c2, _dt, _fmt2v,
        thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
        thrust::raw_pointer_cast(_device_data->_d_mass.data()),
        thrust::raw_pointer_cast(_device_data->_d_fx.data()),
        thrust::raw_pointer_cast(_device_data->_d_fy.data()),
        thrust::raw_pointer_cast(_device_data->_d_fz.data()),
        thrust::raw_pointer_cast(_device_data->_d_vx.data()),
        thrust::raw_pointer_cast(_device_data->_d_vy.data()),
        thrust::raw_pointer_cast(_device_data->_d_vz.data()));
}
void DefaultVelocityController::Update3() {
    op::UpdateVelocityOp3<device::DEVICE_GPU>()
    (
        *(_structure_info_data->_num_atoms), _c3, _dt, _fmt2v,
        thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
        thrust::raw_pointer_cast(_device_data->_d_mass.data()),
        thrust::raw_pointer_cast(_device_data->_d_fx.data()),
        thrust::raw_pointer_cast(_device_data->_d_fy.data()),
        thrust::raw_pointer_cast(_device_data->_d_fz.data()),
        thrust::raw_pointer_cast(_device_data->_d_vx.data()),
        thrust::raw_pointer_cast(_device_data->_d_vy.data()),
        thrust::raw_pointer_cast(_device_data->_d_vz.data()));
}
void DefaultVelocityController::Update4() {
    op::UpdateVelocityOp4<device::DEVICE_GPU>()
    (
        *(_structure_info_data->_num_atoms), _c4, _dt, _fmt2v,
        thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
        thrust::raw_pointer_cast(_device_data->_d_mass.data()),
        thrust::raw_pointer_cast(_device_data->_d_fx.data()),
        thrust::raw_pointer_cast(_device_data->_d_fy.data()),
        thrust::raw_pointer_cast(_device_data->_d_fz.data()),
        thrust::raw_pointer_cast(_device_data->_d_vx.data()),
        thrust::raw_pointer_cast(_device_data->_d_vy.data()),
        thrust::raw_pointer_cast(_device_data->_d_vz.data()));
}