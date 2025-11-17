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
  _structure_info_data = DataManager::getInstance().getMDData()->_structure_info_data;
  _device_data = DataManager::getInstance().getDeviceData();
  _box = DataManager::getInstance().getMDData()->_box;
}

GroupController::~GroupController()
{

}

void GroupController::Init() {

    rbmd::Id num_atoms = *(_structure_info_data->_num_atoms);
    std::vector<rbmd::Id> all_atoms(num_atoms);
    std::iota(all_atoms.begin(), all_atoms.end(), 0);
    _groups["all"] = all_atoms;
}

void GroupController::ComputeVCM(const std::string& group_name, Real3 vcm_out) {

    // if (_groups.find(group_name) == _groups.end()) {
    //     // 在Logger中添加错误处理
    //     Logger::Instance().error("Group {} not found!", group_name);
    //     vcm_out.x = vcm_out.y = vcm_out.z= 0.0;
    //     return;
    // }

   // [mass, mom_vx, mom_vy, mom_vz]
    thrust::device_vector<rbmd::Real> d_vcm_contrib(4, 0.0);
    op::ComputeVCMOp<device::DEVICE_GPU>()(
        *(_structure_info_data->_num_atoms),
        thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
        thrust::raw_pointer_cast(_device_data->_d_mass.data()),
        thrust::raw_pointer_cast(_device_data->_d_vx.data()),
        thrust::raw_pointer_cast(_device_data->_d_vy.data()),
        thrust::raw_pointer_cast(_device_data->_d_vz.data()),
        thrust::raw_pointer_cast(d_vcm_contrib.data()));

     thrust::host_vector<rbmd::Real> h_vcm_contrib = d_vcm_contrib;

  // std::cout <<   "h_vcm_contrib :  "<< h_vcm_contrib[0]  << ", " <<h_vcm_contrib[1]  << ", "
  //   << h_vcm_contrib[2] << ", " <<h_vcm_contrib[2] << std::endl;


    rbmd::Real total_mass = h_vcm_contrib[0];
    if (total_mass > 1e-9) { //
        vcm_out.x = h_vcm_contrib[1] / total_mass;
        vcm_out.y  = h_vcm_contrib[2] / total_mass;
        vcm_out.z = h_vcm_contrib[3] / total_mass;
    } else {
        vcm_out.x = vcm_out.y = vcm_out.z = 0.0;
    }
}

void GroupController::ComputeXCM(const std::string& group_name, Real3 xcm_out) {
   // if (_groups.find(group_name) == _groups.end())
   // { /* error handling */
   //   return;
   // }

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


  // d_xcm_contrib:  [mass, mom_px, mom_py, mom_pz]
  thrust::device_vector<rbmd::Real> d_xcm_contrib(4, 0.0);
  op::ComputeXCMOp<device::DEVICE_GPU>()(
      num_atoms,thrust::raw_pointer_cast(_device_data->_d_mass.data()),
      thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
      thrust::raw_pointer_cast(_device_data->_d_unwarp_px.data()),
      thrust::raw_pointer_cast(_device_data->_d_unwarp_py.data()),
      thrust::raw_pointer_cast(_device_data->_d_unwarp_pz.data()),
      thrust::raw_pointer_cast(d_xcm_contrib.data()));

  thrust::host_vector<rbmd::Real> h_xcm_contrib = d_xcm_contrib;

  rbmd::Real total_mass = h_xcm_contrib[0];
  if (total_mass > 1e-9) {
    xcm_out.x = h_xcm_contrib[1] / total_mass;
    xcm_out.y= h_xcm_contrib[2] / total_mass;
    xcm_out.z = h_xcm_contrib[3] / total_mass;
  } else {
    xcm_out.x= xcm_out.y = xcm_out.z = 0.0;
  }

  // std::cout <<   "h_xcm_contrib :  "<< h_xcm_contrib[0]  << ", " <<h_xcm_contrib[1]  << ", "
  // << h_xcm_contrib[2] << ", " <<h_xcm_contrib[2] << std::endl;
}

void GroupController::ComputeAngMom(const std::string& group_name, Real3 cm, rbmd::Real* angmom_out)
{
  // if (_groups.find(group_name) == _groups.end()) { /* error handling */ return; }

  // 调用GPU Op计算角动量
  //d_angmom_contrib : [Lx, Ly, Lz]
  thrust::device_vector<rbmd::Real> d_angmom_contrib(3, 0.0);
  op::ComputeAngMomOp<device::DEVICE_GPU>()(
      *(_structure_info_data->_num_atoms), cm,
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
      *_box,thrust::raw_pointer_cast(d_angmom_contrib.data()));

  CHECK_RUNTIME(MEMCPY(angmom_out, thrust::raw_pointer_cast(d_angmom_contrib.data()), 3 * sizeof(rbmd::Real), D2H));
  // thrust::host_vector<rbmd::Real> h_angmom_contri = d_angmom_contrib;
  // std::cout <<   "h_angmom_contri :  "<< h_angmom_contri[0]  << ", " <<h_angmom_contri[1]  << ", "<< h_angmom_contri[2]  << std::endl;
}

void GroupController::ComputeInertia(const std::string& group_name, Real3 cm, rbmd::Real (*inertia_out)[3]) {
  // if (_groups.find(group_name) == _groups.end()) { /* error handling */ return; }

  //d_inertia_contrib : [Ixx, Iyy, Izz, Ixy, Ixz, Iyz]
  thrust::device_vector<rbmd::Real> d_inertia_contrib(6, 0.0);
  op::ComputeInertiaOp<device::DEVICE_GPU>()(
       *(_structure_info_data->_num_atoms), cm,
      thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
      thrust::raw_pointer_cast(_device_data->_d_mass.data()),
      thrust::raw_pointer_cast(_device_data->_d_px.data()),
      thrust::raw_pointer_cast(_device_data->_d_py.data()),
      thrust::raw_pointer_cast(_device_data->_d_pz.data()),
      thrust::raw_pointer_cast(_device_data->_d_flagX.data()),
      thrust::raw_pointer_cast(_device_data->_d_flagY.data()),
      thrust::raw_pointer_cast(_device_data->_d_flagZ.data()),
      *_box,thrust::raw_pointer_cast(d_inertia_contrib.data()));

  rbmd::Real h_inertia_contrib[6];
  CHECK_RUNTIME(MEMCPY(h_inertia_contrib, thrust::raw_pointer_cast(d_inertia_contrib.data()), 6 * sizeof(rbmd::Real), D2H));

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

    if (ABS(determinant) > epsilon) {
        // Non-singular matrix: Manually calculate the inverse matrix I^-1
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
      // Use the Jacobi iteration method to solve the eigenvalues (idiag) and eigenvectors (evectors)
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
        //
        Logger::Instance().error("Insufficient Jacobi rotations for omega");
        exit(EXIT_FAILURE); //
      }

      // Extract eigenvectors
      rbmd::Real ex[3] = {evectors[0][0], evectors[1][0], evectors[2][0]};
      rbmd::Real ey[3] = {evectors[0][1], evectors[1][1], evectors[2][1]};
      rbmd::Real ez[3] = {evectors[0][2], evectors[1][2], evectors[2][2]};

      // Check and enforce that the eigenvectors constitutes a right-handed coordinate system
      rbmd::Real cross_product_check[3];
      MathLib::cross3(ex, ey, cross_product_check);
      if (MathLib::dot3(cross_product_check, ez) < 0.0) {
        MathLib::negate3(ez);
      }

      // Check whether the principal inertia (eigenvalue) is close to zero. If so, set it to zero.
      rbmd::Real max_diag;
      max_diag = MAX(idiag[0],idiag[1]);
      max_diag = MAX(max_diag,idiag[2]);

      if (idiag[0] < epsilon * max_diag) idiag[0] = 0.0;
      if (idiag[1] < epsilon * max_diag) idiag[1] = 0.0;
      if (idiag[2] < epsilon * max_diag) idiag[2] = 0.0;

      // Finally, the angular velocity is calculated using the inertia tensor after diagonalization.
      MathLib::angmom_to_omega(angmom, ex, ey, ez, idiag, omega_out);
    }
}