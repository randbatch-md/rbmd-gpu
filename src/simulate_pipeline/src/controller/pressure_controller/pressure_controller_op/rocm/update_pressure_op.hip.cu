#include <hipcub/hipcub.hpp>

#include "rbmd_define.h"
#include "model/box.h"
#include "update_pressure_op.h"

namespace op {

__global__ void X2Lamda(Box* box, const rbmd::Id num_atoms,rbmd::Real* px,
   rbmd::Real* py, rbmd::Real* pz)
{
  unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid1 < num_atoms) {
    Real3 delta;
    delta.data[0] = px[tid1]- box->_coord_min[0];
    delta.data[1] = py[tid1]- box->_coord_min[1];
    delta.data[2] = pz[tid1]- box->_coord_min[2];

    Real3 lamda_position;
    lamda_position.data[0] = box->_length_inv[0] * delta.data[0] + box->_length_inv[5]
    * delta.data[1] + box->_length_inv[4] * delta.data[2];
    lamda_position.data[1] = box->_length_inv[1] * delta.data[1] + box->_length_inv[3]
    * delta.data[2];
    lamda_position.data[2] = box->_length_inv[2] * delta.data[2];

    px[tid1] =  lamda_position.data[0];
    py[tid1] =  lamda_position.data[1];
    pz[tid1] =  lamda_position.data[2];
  }
}

__global__ void Lamda2X(Box* box, const rbmd::Id num_atoms,rbmd::Real* px,
   rbmd::Real* py, rbmd::Real* pz) {
  unsigned int tid1 = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid1 < num_atoms) {
    Real3 position_base;
    position_base.data[0] = px[tid1];
    position_base.data[1] = py[tid1];
    position_base.data[2] = pz[tid1];

    // compute  real position
    Real3 x_position;
    x_position.data[0] = box->_length[0] * position_base.data[0] +
      box->_length[5] * position_base.data[1] +box->_length[4] *
        position_base.data[2] + box->_coord_min[0];

    x_position.data[1] = box->_length[1] * position_base.data[1] +
      box->_length[3] * position_base.data[2] +box->_coord_min[1];

    x_position.data[2] = box->_length[2] * position_base.data[2] + box->_coord_min[2];

    px[tid1]= x_position.data[0];
    py[tid1]= x_position.data[1];
    pz[tid1]= x_position.data[2];

  }
}




void X2LamdaOp<device::DEVICE_GPU>::operator()(
    Box* box, const rbmd::Id num_atoms,rbmd::Real* px,  rbmd::Real* py,
    rbmd::Real* pz ) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

  CHECK_KERNEL(X2Lamda<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>
    (box, num_atoms,px, py, pz));
}

void Lamda2XOp<device::DEVICE_GPU>::operator()(
    Box* box, const rbmd::Id num_atoms,rbmd::Real* px,  rbmd::Real* py,
    rbmd::Real* pz ) {
  unsigned int blocks_per_grid = (num_atoms + BLOCK_SIZE - 1) / BLOCK_SIZE;

  CHECK_KERNEL(Lamda2X<<<blocks_per_grid, BLOCK_SIZE, 0, 0>>>
    (box, num_atoms,px, py, pz));
}

}


