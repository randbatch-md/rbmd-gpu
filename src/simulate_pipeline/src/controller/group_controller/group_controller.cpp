#include "group_controller.h"

#include <thrust/device_ptr.h>

#include <numeric>
#include "math_utils.h" // 包含上面实现的辅助函数
#include "../output/include/Logger.hpp"
#include "device_types.h"
#include "update_temperature_op.h"  //
#include "group_controller_op.h"

GroupController::GroupController() {
    //
  CHECK_RUNTIME(MALLOC(&_d_xcm_contrib, 4 * sizeof(rbmd::Real)));
  CHECK_RUNTIME(MALLOC(&_d_angmom_contrib, 3 * sizeof(rbmd::Real)));
  CHECK_RUNTIME(MALLOC(&_d_inertia_contrib, 6 * sizeof(rbmd::Real)));
}

GroupController::~GroupController() {
  // 在这里添加
  CHECK_RUNTIME(FREE(_d_xcm_contrib));
  CHECK_RUNTIME(FREE(_d_angmom_contrib));
  CHECK_RUNTIME(FREE(_d_inertia_contrib));
}

void GroupController::fetchData() {
    if (!_device_data) {
        _device_data = DataManager::getInstance().getDeviceData();
    }
    if (!_structure_info_data) {
        _structure_info_data = DataManager::getInstance().getMDData()->_structure_info_data;
    }
}

void GroupController::Init() {
    if (_is_initialized) return;

    fetchData();
    
    // 当前只实现一个"all"组
    rbmd::Id num_atoms = *(_structure_info_data->_num_atoms);
    std::vector<rbmd::Id> all_atoms(num_atoms);
    std::iota(all_atoms.begin(), all_atoms.end(), 0);
    _groups["all"] = all_atoms;

    _is_initialized = true;
}

void GroupController::ComputeVCM(const std::string& group_name, Real3 vcm_out) {
    
    if (_groups.find(group_name) == _groups.end()) {
        // 在Logger中添加错误处理
        Logger::Instance().error("Group {} not found!", group_name);
        vcm_out.x = vcm_out.y = vcm_out.z= 0.0;
        return;
    }

    rbmd::Id num_atoms = *(_structure_info_data->_num_atoms);
    _d_vcm_contrib.resize(4,0.0);

    // 调用GPU Op计算总质量和总动量
    op::ComputeVCMOp<device::DEVICE_GPU>()(
        num_atoms,
        thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
        thrust::raw_pointer_cast(_device_data->_d_mass.data()),
        thrust::raw_pointer_cast(_device_data->_d_vx.data()),
        thrust::raw_pointer_cast(_device_data->_d_vy.data()),
        thrust::raw_pointer_cast(_device_data->_d_vz.data()),
        thrust::raw_pointer_cast(_d_vcm_contrib.data()));

     thrust::host_vector<rbmd::Real> h_vcm_contrib;
     h_vcm_contrib.resize(4,0.0);
     h_vcm_contrib = _d_vcm_contrib;

  // std::cout <<   "h_vcm_contrib :  "<< h_vcm_contrib[0]  << ", " <<h_vcm_contrib[1]  << ", "
  //   << h_vcm_contrib[2] << ", " <<h_vcm_contrib[2] << std::endl;


    rbmd::Real total_mass = h_vcm_contrib[0];
    if (total_mass > 1e-9) { // 避免除以零
        vcm_out.x = h_vcm_contrib[1] / total_mass;
        vcm_out.y  = h_vcm_contrib[2] / total_mass;
        vcm_out.z = h_vcm_contrib[3] / total_mass;
    } else {
        vcm_out.x = vcm_out.y = vcm_out.z = 0.0;
    }
}

void GroupController::ComputeXCM(const std::string& group_name, Real3 xcm_out) {
  if (!_is_initialized) Init();
  if (_groups.find(group_name) == _groups.end())
  { /* error handling */
    return;
  }

  //Unwarp  position
  rbmd::Id num_atoms = *(_structure_info_data->_num_atoms);
  op::UnwarpPositionOp<device::DEVICE_GPU>()(num_atoms,*_box,
    thrust::raw_pointer_cast(_device_data->_d_px.data()),
    thrust::raw_pointer_cast(_device_data->_d_py.data()),
    thrust::raw_pointer_cast(_device_data->_d_pz.data()),
    thrust::raw_pointer_cast(_device_data->_d_flagX.data()),
    thrust::raw_pointer_cast(_device_data->_d_flagY.data()),
    thrust::raw_pointer_cast(_device_data->_d_flagZ.data()),
    thrust::raw_pointer_cast(_device_data->_d_unwarp_px.data()),
    thrust::raw_pointer_cast(_device_data->_d_unwarp_py.data()),
    thrust::raw_pointer_cast(_device_data->_d_unwarp_pz.data()));


  // 调用GPU Op计算总质量和质量矩
  CHECK_RUNTIME(MEMSET(_d_xcm_contrib, 0, 4 * sizeof(rbmd::Real)));
  op::ComputeXCMOp<device::DEVICE_GPU>()(
      num_atoms,thrust::raw_pointer_cast(_device_data->_d_mass.data()),
      thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
      thrust::raw_pointer_cast(_device_data->_d_unwarp_px.data()),
      thrust::raw_pointer_cast(_device_data->_d_unwarp_py.data()),
      thrust::raw_pointer_cast(_device_data->_d_unwarp_pz.data()),
      _d_xcm_contrib);

  rbmd::Real h_xcm_contrib[4];
  CHECK_RUNTIME(MEMCPY(h_xcm_contrib, _d_xcm_contrib, 4 * sizeof(rbmd::Real), D2H));

  rbmd::Real total_mass = h_xcm_contrib[0];
  if (total_mass > 1e-9) {
    xcm_out.x = h_xcm_contrib[1] / total_mass;
    xcm_out.y= h_xcm_contrib[2] / total_mass;
    xcm_out.z = h_xcm_contrib[3] / total_mass;
  } else {
    xcm_out.x= xcm_out.y = xcm_out.z = 0.0;
  }
}

void GroupController::ComputeAngMom(const std::string& group_name, Real3 cm, rbmd::Real* angmom_out)
{
  if (!_is_initialized) Init();
  if (_groups.find(group_name) == _groups.end()) { /* error handling */ return; }

  rbmd::Id num_atoms = *(_structure_info_data->_num_atoms);
  CHECK_RUNTIME(MEMSET(_d_angmom_contrib, 0, 3 * sizeof(rbmd::Real)));

  // 调用GPU Op计算角动量
  op::ComputeAngMomOp<device::DEVICE_GPU>()(
      num_atoms, cm,
      thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
      thrust::raw_pointer_cast(_device_data->_d_mass.data()),
      thrust::raw_pointer_cast(_device_data->_d_px.data()),
      thrust::raw_pointer_cast(_device_data->_d_py.data()),
      thrust::raw_pointer_cast(_device_data->_d_pz.data()),
      thrust::raw_pointer_cast(_device_data->_d_vx.data()),
      thrust::raw_pointer_cast(_device_data->_d_vy.data()),
      thrust::raw_pointer_cast(_device_data->_d_vz.data()),
      thrust::raw_pointer_cast(_device_data->_d_flagX.data()),
      thrust::raw_pointer_cast(_device_data->_d_flagY.data()),
      thrust::raw_pointer_cast(_device_data->_d_flagZ.data()),
      *_box,_d_angmom_contrib);

  CHECK_RUNTIME(MEMCPY(angmom_out, _d_angmom_contrib, 3 * sizeof(rbmd::Real), D2H));
}

void GroupController::ComputeInertia(const std::string& group_name, Real3 cm, rbmd::Real (*inertia_out)[3]) {
  if (!_is_initialized) Init();
  if (_groups.find(group_name) == _groups.end()) { /* error handling */ return; }

  rbmd::Id num_atoms = *(_structure_info_data->_num_atoms);
  CHECK_RUNTIME(MEMSET(_d_inertia_contrib, 0, 6 * sizeof(rbmd::Real)));

  op::ComputeInertiaOp<device::DEVICE_GPU>()(
      num_atoms, cm,
      thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
      thrust::raw_pointer_cast(_device_data->_d_mass.data()),
      thrust::raw_pointer_cast(_device_data->_d_px.data()),
      thrust::raw_pointer_cast(_device_data->_d_py.data()),
      thrust::raw_pointer_cast(_device_data->_d_pz.data()),
      thrust::raw_pointer_cast(_device_data->_d_flagX.data()),
      thrust::raw_pointer_cast(_device_data->_d_flagY.data()),
      thrust::raw_pointer_cast(_device_data->_d_flagZ.data()),
      *_box,_d_inertia_contrib);

  rbmd::Real h_inertia_contrib[6];
  CHECK_RUNTIME(MEMCPY(h_inertia_contrib, _d_inertia_contrib, 6 * sizeof(rbmd::Real), D2H));

  inertia_out[0][0] = h_inertia_contrib[0];
  inertia_out[1][1] = h_inertia_contrib[1];
  inertia_out[2][2] = h_inertia_contrib[2];
  inertia_out[0][1] = inertia_out[1][0] = h_inertia_contrib[3];
  inertia_out[0][2] = inertia_out[2][0] = h_inertia_contrib[4];
  inertia_out[1][2] = inertia_out[2][1] = h_inertia_contrib[5];
}

void GroupController::ComputeOmega(const rbmd::Real* angmom, const rbmd::Real (*inertia)[3], rbmd::Real* omega_out) {
    const rbmd::Real epsilon  = 1.0e-6;

    const rbmd::Real determinant = inertia[0][0] * inertia[1][1] * inertia[2][2]
                                 + inertia[0][1] * inertia[1][2] * inertia[2][0]
                                 + inertia[0][2] * inertia[1][0] * inertia[2][1]
                                 - inertia[0][0] * inertia[1][2] * inertia[2][1]
                                 - inertia[0][1] * inertia[1][0] * inertia[2][2]
                                 - inertia[2][0] * inertia[1][1] * inertia[0][2];

    if (std::fabs(determinant) > epsilon) {
        // 非奇异矩阵: 手动计算逆矩阵 I^-1
        rbmd::Real inverse[3][3];
        inverse[0][0] = inertia[1][1] * inertia[2][2] - inertia[1][2] * inertia[2][1];
        inverse[0][1] = -(inertia[0][1] * inertia[2][2] - inertia[0][2] * inertia[2][1]);
        inverse[0][2] = inertia[0][1] * inertia[1][2] - inertia[0][2] * inertia[1][1];
        inverse[1][0] = -(inertia[1][0] * inertia[2][2] - inertia[1][2] * inertia[2][0]);
        inverse[1][1] = inertia[0][0] * inertia[2][2] - inertia[0][2] * inertia[2][0];
        inverse[1][2] = -(inertia[0][0] * inertia[1][2] - inertia[0][2] * inertia[1][0]);
        inverse[2][0] = inertia[1][0] * inertia[2][1] - inertia[1][1] * inertia[2][0];
        inverse[2][1] = -(inertia[0][0] * inertia[2][1] - inertia[0][1] * inertia[2][0]);
        inverse[2][2] = inertia[0][0] * inertia[1][1] - inertia[0][1] * inertia[1][0];

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                inverse[i][j] /= determinant;
            }
        }

        // omega = I^-1 * L
        omega_out[0] = inverse[0][0] * angmom[0] + inverse[0][1] * angmom[1] + inverse[0][2] * angmom[2];
        omega_out[1] = inverse[1][0] * angmom[0] + inverse[1][1] * angmom[1] + inverse[1][2] * angmom[2];
        omega_out[2] = inverse[2][0] * angmom[0] + inverse[2][1] * angmom[1] + inverse[2][2] * angmom[2];

    } else {
      // 使用雅可比迭代法求解特征值（idiag）和特征向量（evectors）
      rbmd::Real idiag[3];
      rbmd::Real evectors[3][3];
      rbmd::Real inertia_copy[3][3];

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          inertia_copy[i][j] = inertia[i][j];
        }
      }

      int ierror = MathLib::jacobi3(inertia_copy, idiag, evectors);
      if (ierror) {
        // 如果未收敛，抛出错误或警告
        Logger::Instance().error("Insufficient Jacobi rotations for omega");
        exit(EXIT_FAILURE); //
      }

      // 提取特征向量
      rbmd::Real ex[3] = {evectors[0][0], evectors[1][0], evectors[2][0]};
      rbmd::Real ey[3] = {evectors[0][1], evectors[1][1], evectors[2][1]};
      rbmd::Real ez[3] = {evectors[0][2], evectors[1][2], evectors[2][2]};

      // 检查并强制特征向量构成右手坐标系
      rbmd::Real cross_product_check[3];
      MathLib::cross3(ex, ey, cross_product_check);
      if (MathLib::dot3(cross_product_check, ez) < 0.0) {
        MathLib::negate3(ez);
      }

      // 检查主惯量（特征值）是否接近于零，如果是，将其置零
      rbmd::Real max_diag;
      max_diag = MAX(idiag[0],idiag[1]);
      max_diag = MAX(max_diag,idiag[2]);

      if (idiag[0] < epsilon * max_diag) idiag[0] = 0.0;
      if (idiag[1] < epsilon * max_diag) idiag[1] = 0.0;
      if (idiag[2] < epsilon * max_diag) idiag[2] = 0.0;

      // 最后，使用对角化后的惯量张量计算角速度
      MathLib::angmom_to_omega(angmom, ex, ey, ez, idiag, omega_out);
    }
}