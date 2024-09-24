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

void LJForce::Init() { _num_atoms = *(_structure_info_data->_num_atoms);}

void LJForce::Execute() {

  rbmd::Real cut_off = 5.0;
  //
  extern int test_current_step;
  list = full_list_builder->Build();
  //list->print("./neighbor_list_step_" + std::to_string(test_current_step)+"_.csv");

  rbmd::Real h_total_evdwl = 0.0;

  CHECK_RUNTIME(MEMSET(_d_total_evdwl,0,sizeof(rbmd::Real)));
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
              _d_total_evdwl);


  CHECK_RUNTIME(MEMCPY(&h_total_evdwl,_d_total_evdwl , sizeof(rbmd::Real), D2H));

  // 打印累加后的总能量
  rbmd::Real  ave_evdwl = h_total_evdwl/_num_atoms;
  std::cout << "test_current_step:" << test_current_step <<  " " << "average_vdwl_energy:" << ave_evdwl << std::endl;



  std::cout << "out of force execute" << std::endl;


  //out
  std::ofstream outfile("ave_evdwl.txt", std::ios::app);
  outfile << test_current_step << " " << ave_evdwl << std::endl;
  outfile.close();
             
}

