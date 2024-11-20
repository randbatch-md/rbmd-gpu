#include <hipcub/hipcub.hpp>

#include "rbmd_define.h"
#include "model/box.h"
#include "update_pressure_op.h"

namespace op {

__global__ void X2Lamda(Box box, const rbmd::Id num_atoms,rbmd::Real* px,
   rbmd::Real* py, rbmd::Real* pz)
{
  unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid1 < num_atoms) {
    Real3 delta;
    REAL_DATA(delta)[0] = px[tid1]- box._coord_min[0];
    REAL_DATA(delta)[1] = py[tid1]- box._coord_min[1];
    REAL_DATA(delta)[2] = pz[tid1]- box._coord_min[2];

    Real3 lamda_position;
    REAL_DATA(lamda_position)[0] = box._length_inv[0] * REAL_DATA(delta)[0] + box._length_inv[5]
    * REAL_DATA(delta)[1] + box._length_inv[4] * REAL_DATA(delta)[2];
    REAL_DATA(lamda_position)[1] = box._length_inv[1] * REAL_DATA(delta)[1] + box._length_inv[3]
    * REAL_DATA(delta)[2];
    REAL_DATA(lamda_position)[2] = box._length_inv[2] * REAL_DATA(delta)[2];

    px[tid1] =  REAL_DATA(lamda_position)[0];
    py[tid1] =  REAL_DATA(lamda_position)[1];
    pz[tid1] =  REAL_DATA(lamda_position)[2];
  }
}

__global__ void Lamda2X(Box box, const rbmd::Id num_atoms,rbmd::Real* px,
   rbmd::Real* py, rbmd::Real* pz) {
  unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid1 < num_atoms) {
    Real3 position_base;
    REAL_DATA(position_base)[0] = px[tid1];
    REAL_DATA(position_base)[1] = py[tid1];
    REAL_DATA(position_base)[2] = pz[tid1];

    // compute  real position
    Real3 x_position;
    REAL_DATA(x_position)[0] = box._length[0] * REAL_DATA(position_base)[0] +
      box._length[5] * REAL_DATA(position_base)[1] +box._length[4] *
        REAL_DATA(position_base)[2] + box._coord_min[0];

    REAL_DATA(x_position)[1] = box._length[1] * REAL_DATA(position_base)[1] +
      box._length[3] * REAL_DATA(position_base)[2] +box._coord_min[1];

    REAL_DATA(x_position)[2] = box._length[2] * REAL_DATA(position_base)[2] + box._coord_min[2];

    px[tid1]= REAL_DATA(x_position)[0];
    py[tid1]= REAL_DATA(x_position)[1];
    pz[tid1]= REAL_DATA(x_position)[2];

  }
}

__global__ void UpdataVelocityRescalePressure(const rbmd::Id num_atoms,
                                      const Real3 factor,
                                      rbmd::Real* vx, rbmd::Real* vy,
                                      rbmd::Real* vz) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  if (tid < num_atoms) {
    vx[tid] *= factor.x;
    vy[tid] *= factor.y;
    vz[tid] *= factor.z;
  }
}


void X2LamdaOp<device::DEVICE_GPU>::operator()(
   Box box, const rbmd::Id num_atoms,rbmd::Real* px,  rbmd::Real* py,
    rbmd::Real* pz ) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

  CHECK_KERNEL(X2Lamda<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>
    (box, num_atoms,px, py, pz));
}

void Lamda2XOp<device::DEVICE_GPU>::operator()(
   Box box, const rbmd::Id num_atoms,rbmd::Real* px,  rbmd::Real* py,
    rbmd::Real* pz ) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

  CHECK_KERNEL(Lamda2X<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>
    (box, num_atoms,px, py, pz));
}

void UpdataVelocityRescalePressureOp<device::DEVICE_GPU>::operator()(
    const rbmd::Id num_atoms, const Real3 factor, rbmd::Real* vx,
    rbmd::Real* vy, rbmd::Real* vz) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;
  CHECK_KERNEL(UpdataVelocityRescalePressure<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>(
      num_atoms, factor, vx, vy, vz));

}
}


