#include "ljforce.h"

#include <thrust/device_ptr.h>

#include "../../common/device_types.h"
#include "../../common/types.h"
#include "ljforce_op/ljforce_op.h"
// #include <hipcub/hipcub.hpp>
// #include <hipcub/backend/rocprim/block/block_reduce.hpp>

LJForce::LJForce() {
  full_list_builder = std::make_shared<FullNeighborListBuilder>();
  list = full_list_builder->Build();
};

void LJForce::Init() { _num_atoms = *(_structure_info_data->_num_atoms); }

void LJForce::Execute() {
  LJForce::Init();
  rbmd::Real cut_off = 5.0;
//   std::cout << "  _h_total_max_neighbor_num   "
// << list->_h_total_max_neighbor_num << std::endl;
//   std::cout << " _d_neighbor_num 0  " << list->_d_neighbor_num[0] << std::endl;
//   std::cout << " _d_neighbor_num 20  " << list->_d_neighbor_num[20] << std::endl;
//   std::cout << "0 邻居起始索引：" << list->_start_idx[0] << "  0 邻居结束索引："
//         << list->_end_idx[0] << std::endl;
//   std::cout << "1 邻居起始索引：" << list->_start_idx[1] << "  1 邻居结束索引："
//           << list->_end_idx[1] << std::endl;
//   std::cout << "1 的第一个邻居：" << list->_d_neighbors[list->_start_idx[1]]
//             << std::endl;
//   std::cout << "1 的最后一个邻居：" << list->_d_neighbors[list->_end_idx[1]]
//           << std::endl;
//   std::cout << "===================       原子idx20的邻居索引：         =====================" << std::endl;
//   for (int i = list->_start_idx[20]; i < list->_end_idx[20]; ++i)
// {
//   std::cout << "  " << list->_d_neighbors[i] <<std::endl;
// }
  op::LJForceOp<device::DEVICE_GPU> lj_force_op;
  lj_force_op(_device_data->_d_box, cut_off, _num_atoms,
              thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
              thrust::raw_pointer_cast(_device_data->_d_molecular_id.data()),
              thrust::raw_pointer_cast(_device_data->_d_eps.data()),
              thrust::raw_pointer_cast(_device_data->_d_sigma.data()),
              thrust::raw_pointer_cast(list->_start_idx.data()),
              thrust::raw_pointer_cast(list->_end_idx.data()),
              thrust::raw_pointer_cast(list->_d_neighbors.data()),
              thrust::raw_pointer_cast(_device_data->_d_px.data()),
              thrust::raw_pointer_cast(_device_data->_d_py.data()),
              thrust::raw_pointer_cast(_device_data->_d_pz.data()),
              thrust::raw_pointer_cast(_device_data->_d_fx.data()),
              thrust::raw_pointer_cast(_device_data->_d_fy.data()),
              thrust::raw_pointer_cast(_device_data->_d_fz.data()),
              thrust::raw_pointer_cast(_device_data->_d_evdwl.data()));

  std::cout << "out of force execute" << std::endl;
  // // ����������
  // rbmd::Real* d_total_evdwl;
  // hipMalloc(&d_total_evdwl, sizeof(rbmd::Real));
  // op::DeviceReduceSum(thrust::raw_pointer_cast(_device_data->_d_evdwl.data()),
  // d_total_evdwl, _num_atoms);

  //// ���豸���ƽ��������
  // rbmd::Real total_evdwl;
  // hipMemcpy(&total_evdwl, d_total_evdwl, sizeof(rbmd::Real),
  // hipMemcpyDeviceToHost);

  // std::cout << "Total evdwl: " << total_evdwl << std::endl;

    //hipFree(d_total_evdwl);

             
}

