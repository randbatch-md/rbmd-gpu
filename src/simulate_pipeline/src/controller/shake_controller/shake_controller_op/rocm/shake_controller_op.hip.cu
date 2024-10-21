#include <hip/hip_runtime.h>

#include "rbmd_define.h"
#include "shake_controller_op.h"
#include <math.h>
namespace op {
#define THREADS_PER_BLOCK 256

__device__ Real3 MinDistanceVec(const rbmd::Real& share_px_1,
                                const rbmd::Real& share_py_1, 
                                const rbmd::Real& share_pz_1, 
                                const rbmd::Real& share_px_2,
                                const rbmd::Real& share_py_2,
                                const rbmd::Real& share_pz_2,
                                Box* box)
{
    rbmd::Id periodicX = 1;
    rbmd::Id periodicY = 1;
    rbmd::Id periodicZ = 1;

    Real3 vec;
    vec.x = share_px_1 - share_px_2;
    vec.y = share_py_1 - share_py_2;
    vec.z = share_pz_1 - share_pz_2;
    
    // X
    if (periodicX)
    {
        if (ABS(vec.x) > box->_length[0] * 0.5)
        {
            vec.x -= (vec.x > 0 ? box->_length[0] : -(box->_length[0]));
        }
    }

    // Y
    if (periodicY)
    {
        if (ABS(vec.y) > box->_length[1] * 0.5)
        {
            vec.y -= (vec.y > 0 ? box->_length[1] : -(box->_length[1]));
        }
    }

    // Z
    if (periodicZ)
    {
        if (ABS(vec.z) > box->_length[2] * 0.5)
        {
            vec.z -= (vec.z > 0 ? box->_length[2] : -(box->_length[2]));
        }
    }

    return vec;
}

__global__ void ShakeA(const rbmd::Id num_angle,
                       const rbmd::Real dt,
                       const rbmd::Real fmt2v,
                       Box* box,
                       const rbmd::Real min_x, const rbmd::Real min_y,
                       const rbmd::Real min_z, const rbmd::Real max_x,
                       const rbmd::Real max_y, const rbmd::Real max_z,
                       const rbmd::Real* mass,
                       const rbmd::Id* atoms_type,
                       const Id3* angle_id_vec,
                       rbmd::Real* shake_px,
                       rbmd::Real* shake_py,
                       rbmd::Real* shake_pz,
                       rbmd::Real* shake_vx,
                       rbmd::Real* shake_vy,
                       rbmd::Real* shake_vz,
                       const rbmd::Real* fx,
                       const rbmd::Real* fy,
                       const rbmd::Real* fz,
                       rbmd::Id* flag_px,
                       rbmd::Id* flag_py,
                       rbmd::Id* flag_pz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_angle) {
    rbmd::Id id_0 = angle_id_vec[tid].x;
    rbmd::Id id_1 = angle_id_vec[tid].y;
    rbmd::Id id_2 = angle_id_vec[tid].z;

    Real3 shake_position_0, shake_position_1, shake_position_2;
    shake_position_0.x = shake_px[id_0] + shake_vx[id_0] * dt + 0.5 * dt * dt * fx[id_0] / mass[atoms_type[id_0]] * fmt2v;
    shake_position_0.y = shake_py[id_0] + shake_vy[id_0] * dt + 0.5 * dt * dt * fy[id_0] / mass[atoms_type[id_0]] * fmt2v;
    shake_position_0.z = shake_pz[id_0] + shake_vz[id_0] * dt + 0.5 * dt * dt * fz[id_0] / mass[atoms_type[id_0]] * fmt2v;

    shake_position_1.x = shake_px[id_1] + shake_vx[id_1] * dt + 0.5 * dt * dt * fx[id_1] / mass[atoms_type[id_1]] * fmt2v;
    shake_position_1.y = shake_py[id_1] + shake_vy[id_1] * dt + 0.5 * dt * dt * fy[id_1] / mass[atoms_type[id_1]] * fmt2v;
    shake_position_1.z = shake_pz[id_1] + shake_vz[id_1] * dt + 0.5 * dt * dt * fz[id_1] / mass[atoms_type[id_1]] * fmt2v;

    shake_position_2.x = shake_px[id_2] + shake_vx[id_2] * dt + 0.5 * dt * dt * fx[id_2] / mass[atoms_type[id_2]] * fmt2v;
    shake_position_2.y = shake_py[id_2] + shake_vy[id_2] * dt + 0.5 * dt * dt * fy[id_2] / mass[atoms_type[id_2]] * fmt2v;
    shake_position_2.z = shake_pz[id_2] + shake_vz[id_2] * dt + 0.5 * dt * dt * fz[id_2] / mass[atoms_type[id_2]] * fmt2v;

    rbmd::Real bond1 = 1.0;
    rbmd::Real bond2 = 1.0;
    rbmd::Real bond12 = SQRT(bond1 * bond1 + bond2 * bond2 - 2.0 * bond1 * bond2 * COS((109.4700 / 180.0) * M_PIf));

    
    // minimum image
    Real3 r01 = MinDistanceVec(shake_px[id_1], shake_py[id_1], shake_pz[id_1], shake_px[id_0], shake_py[id_0], shake_pz[id_0], box);
    Real3 r12 = MinDistanceVec(shake_px[id_2], shake_py[id_2], shake_pz[id_2], shake_px[id_1], shake_py[id_1], shake_pz[id_1], box);
    Real3 r20 = MinDistanceVec(shake_px[id_0], shake_py[id_0], shake_pz[id_0], shake_px[id_2], shake_py[id_2], shake_pz[id_2], box);


    // s01,s02,s12 = distance vec after unconstrained update, with PBC
    Real3 s10 = MinDistanceVec(shake_position_0.x, shake_position_0.y, shake_position_0.z, shake_position_1.x, shake_position_1.y, shake_position_1.z, box);
    Real3 s21 = MinDistanceVec(shake_position_1.x, shake_position_1.y, shake_position_1.z, shake_position_2.x, shake_position_2.y, shake_position_2.z, box);
    Real3 s02 = MinDistanceVec(shake_position_2.x, shake_position_2.y, shake_position_2.z, shake_position_0.x, shake_position_0.y, shake_position_0.z, box);

  }
}

void ShakeAOp<device::DEVICE_GPU>::operator()(const rbmd::Id num_angle,
                                              const rbmd::Real dt,
                                              const rbmd::Real fmt2v,
                                              Box* box,
                                              const rbmd::Real min_x, const rbmd::Real min_y,
                                              const rbmd::Real min_z, const rbmd::Real max_x,
                                              const rbmd::Real max_y, const rbmd::Real max_z,
                                              const rbmd::Real* mass,
                                              const rbmd::Id* atoms_type,
                                              const Id3* angle_id_vec,
                                              rbmd::Real* shake_px,
                                              rbmd::Real* shake_py,
                                              rbmd::Real* shake_pz,
                                              rbmd::Real* shake_vx,
                                              rbmd::Real* shake_vy,
                                              rbmd::Real* shake_vz,
                                              const rbmd::Real* fx,
                                              const rbmd::Real* fy,
                                              const rbmd::Real* fz,
                                              rbmd::Id* flag_px,
                                              rbmd::Id* flag_py,
                                              rbmd::Id* flag_pz) 
{
  unsigned int blocks_per_grid = (num_angle + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(ShakeA<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_angle, dt, fmt2v, box, min_x, min_y, min_z, max_x, max_y, max_z, mass, atoms_type, angle_id_vec,
      shake_px, shake_py, shake_pz, shake_vx, shake_vy, shake_vz, fx, fy, fz, flag_px, flag_py, flag_pz));
}
}  // namespace op
