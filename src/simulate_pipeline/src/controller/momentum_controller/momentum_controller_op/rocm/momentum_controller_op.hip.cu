#include "../data_manager/include/model/box.h"
#include "momentum_controller_op.h"
#include "rbmd_define.h"

namespace op {
#define THREADS_PER_BLOCK 256

__device__ inline void unwarp_mom(const rbmd::Real* pos, const int* image,
  Box box, rbmd::Real* unwrap_pos)
{

  unwrap_pos[0] = pos[0] + image[0] * box._length[0];
  unwrap_pos[1] = pos[1] + image[1] * box._length[1];
  unwrap_pos[2] = pos[2] + image[2] * box._length[2];
}

/**
 * @brief GPU内核：计算一个原子组的总动能
 */
__global__ void compute_ke_kernel(const rbmd::Id num_atoms,
                                  const rbmd::Id* atoms_type,
                                  const rbmd::Real* mass,
                                  const rbmd::Real* vx,
                                  const rbmd::Real* vy,
                                  const rbmd::Real* vz,
                                  rbmd::Real* ke_contrib) {
  __shared__ typename BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>::TempStorage temp_storage_ke;

  int i = blockIdx.x * blockDim.x + threadIdx.x;
  rbmd::Real local_ke = 0.0;
  if (i < num_atoms) {
    rbmd::Real mass_i = mass[atoms_type[i]];
    rbmd::Real vx_i = vx[i];
    rbmd::Real vy_i = vy[i];
    rbmd::Real vz_i = vz[i];

    rbmd::Real local_ke = mass_i * (vx_i * vx_i + vy_i * vy_i + vz_i * vz_i);
  }

  rbmd::Real block_ke = BLOCKREDUCE<rbmd::Real, BLOCK_SIZE>
   (temp_storage_ke).Sum(local_ke);

  // 线程 0 将块内结果原子性地累加到全局变量
  if (threadIdx.x == 0) {
    atomicAdd(&ke_contrib[0], block_ke);  // 总质量
  }

}


/**
 * @brief GPU内核：通过减去质心速度来移除线动量
 */
__global__ void zero_linear_momentum_kernel(const rbmd::Id num_atoms,
                                            Real3 vcm,
                                            bool x_flag, bool y_flag, bool z_flag,
                                            rbmd::Real* vx,
                                            rbmd::Real* vy,
                                            rbmd::Real* vz) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= num_atoms) return;

    if (x_flag) vx[i] -= vcm.x;
    if (y_flag) vy[i] -= vcm.y;
    if (z_flag) vz[i] -= vcm.z;
}



/**
 * @brief GPU内核：通过减去旋转速度分量来移除角动量
 */
__global__ void zero_angular_momentum_kernel(const rbmd::Id num_atoms,
                                             Real3 xcm,
                                             const rbmd::Real* omega,
                                             const rbmd::Real* px, const rbmd::Real* py, const rbmd::Real* pz,
                                             const int* d_image_x, const int* d_image_y, const int* d_image_z,
                                             Box box,rbmd::Real* d_vx, rbmd::Real* d_vy, rbmd::Real* d_vz)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < num_atoms)
    {
      rbmd::Real pos_wrap[3] = {px[i], py[i], pz[i]};
      int image[3] = {d_image_x[i], d_image_y[i], d_image_z[i]};
      rbmd::Real pos_unwrap[3];
      unwarp_mom(pos_wrap, image, box, pos_unwrap);

      rbmd::Real dx = pos_unwrap[0] - xcm.x;
      rbmd::Real dy = pos_unwrap[1] - xcm.y;
      rbmd::Real dz = pos_unwrap[2] - xcm.z;

      // 计算旋转速度 v_rot = omega X r
      rbmd::Real v_rot_x = omega[1] * dz - omega[2] * dy;
      rbmd::Real v_rot_y = omega[2] * dx - omega[0] * dz;
      rbmd::Real v_rot_z = omega[0] * dy - omega[1] * dx;

      // v_new = v_old - v_rot
      d_vx[i] -= v_rot_x;
      d_vy[i] -= v_rot_y;
      d_vz[i] -= v_rot_z;
    }
}

void ComputeKineticEnergyOp<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms,rbmd::Id* atoms_type, const rbmd::Real* mass,
    const rbmd::Real* vx, const rbmd::Real* vy, const rbmd::Real* vz,
    rbmd::Real* d_ke_contrib) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(compute_ke_kernel<<<blocks_per_grid, BLOCK_SIZE,0,0>>>(num_atoms, atoms_type,
    mass, vx, vy, vz, d_ke_contrib));
}
  void ZeroLinearMomentumOp<device::DEVICE_GPU>::operator()(
    rbmd::Id num_atoms, Real3 vcm,
    bool x_flag, bool y_flag, bool z_flag,
    rbmd::Real* vx, rbmd::Real* vy, rbmd::Real* vz) {

  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(zero_linear_momentum_kernel<<<blocks_per_grid, BLOCK_SIZE,0,0>>>(num_atoms,
    vcm, x_flag, y_flag, z_flag, vx, vy, vz));
}



void ZeroAngularMomentumOp<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms, Real3 xcm, const rbmd::Real* omega,
    const rbmd::Real* px, const rbmd::Real* py, const rbmd::Real* pz,
    const int* d_image_x, const int* d_image_y, const int* d_image_z,
    Box box,rbmd::Real* d_vx, rbmd::Real* d_vy, rbmd::Real* d_vz)
{

  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
   CHECK_KERNEL(zero_angular_momentum_kernel<<<blocks_per_grid, BLOCK_SIZE,0,0>>>(num_atoms, xcm, omega,
    px, py, pz,d_image_x, d_image_y, d_image_z,box, d_vx, d_vy, d_vz));
}

}