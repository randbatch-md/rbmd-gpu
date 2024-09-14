#include "ljforce.h"

#include <thrust/device_ptr.h>

#include "../../common/device_types.h"
#include "../../common/types.h"
#include "ljforce_op/ljforce_op.h"
// #include <hipcub/hipcub.hpp>
// #include <hipcub/backend/rocprim/block/block_reduce.hpp>

LJForce::LJForce()
{
  full_list_builder = std::make_shared<FullNeighborListBuilder>();
};

void LJForce::Init() { _num_atoms = *(_structure_info_data->_num_atoms); }

void LJForce::Execute() {
  LJForce::Init();
  rbmd::Real cut_off = 5.0;
  //

  list = full_list_builder->Build();
  list->print("./step1neighbor.csv");
  //
    rbmd::Real h_total_evdwl = 0.0;
  rbmd::Real* d_total_evdwl;

  CHECK_RUNTIME(MALLOC(&d_total_evdwl, sizeof(rbmd::Real)));
  CHECK_RUNTIME(MEMCPY(d_total_evdwl,&h_total_evdwl , sizeof(rbmd::Real), H2D));

  //compute LJForce
  op::LJForceOp<device::DEVICE_GPU> lj_force_op;
  lj_force_op(_device_data->_d_box, cut_off, _num_atoms,
              thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
              thrust::raw_pointer_cast(_device_data->_d_molecular_id.data()),
              thrust::raw_pointer_cast(_device_data->_d_sigma.data()),
              thrust::raw_pointer_cast(_device_data->_d_eps.data()),
              thrust::raw_pointer_cast(list->_start_idx.data()),
              thrust::raw_pointer_cast(list->_end_idx.data()),
              thrust::raw_pointer_cast(list->_d_neighbors.data()),
              thrust::raw_pointer_cast(_device_data->_d_px.data()),
              thrust::raw_pointer_cast(_device_data->_d_py.data()),
              thrust::raw_pointer_cast(_device_data->_d_pz.data()),
              thrust::raw_pointer_cast(_device_data->_d_fx.data()),
              thrust::raw_pointer_cast(_device_data->_d_fy.data()),
              thrust::raw_pointer_cast(_device_data->_d_fz.data()),
              thrust::raw_pointer_cast(_device_data->_d_evdwl.data()),
              d_total_evdwl);


  CHECK_RUNTIME(MEMCPY(&h_total_evdwl,d_total_evdwl , sizeof(rbmd::Real), D2H));

  // 打印累加后的总能量
  rbmd::Real  ave_evdwl = h_total_evdwl/_num_atoms;
  std::cout << "Total evdwl: " << h_total_evdwl  << ","<<  "ave_evdwl: "  << ave_evdwl << std::endl;

  // 释放设备端分配的内存
  CHECK_RUNTIME(FREE(d_total_evdwl));

     //std::cout <<"id0" << _device_data->_d_atoms_id[0] << std::endl;
     //std::cout << "p_atoms_id[0] "  <<  _device_data->_d_vx[0]<< " "
     //    << _device_data->_d_vy[0] << " "
     //    << _device_data->_d_vz[0]
     //    << std::endl;

  std::cout << "out of force execute" << std::endl;
             
}

