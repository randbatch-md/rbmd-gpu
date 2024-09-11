#include "ljforce.h"
#include "ljforce_op/ljforce_op.h"
#include <thrust/device_ptr.h>
#include "../../common/device_types.h"
#include "../../common/types.h"
//#include <hipcub/hipcub.hpp> 
//#include <hipcub/backend/rocprim/block/block_reduce.hpp>

LJForce::LJForce(){
    full_list_builder = std::make_shared<FullNeighborListBuilder>();
    list = full_list_builder->Build();
};

void LJForce::Init()
{
    _num_atoms = _structure_info_data->_num_atoms;
}

void LJForce::Execute()
{
    LJForce::Init();
    rbmd::Real cut_off = 5.0;

	op::LJForceOp<device::DEVICE_GPU> lj_force_op;

    lj_force_op(_device_data->_d_box,
             cut_off,
             _num_atoms,
             thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
             thrust::raw_pointer_cast(_device_data->_d_molecular_id.data()),
             thrust::raw_pointer_cast(_device_data->_d_eps.data()),
             thrust::raw_pointer_cast(_device_data->_d_sigma.data()),
             thrust::raw_pointer_cast(list->_start_idx.data()),
             thrust::raw_pointer_cast(list->_end_idx.data()),
             thrust::raw_pointer_cast(list->_d_neighbor_num.data()),
             thrust::raw_pointer_cast(_device_data->_d_px.data()),
             thrust::raw_pointer_cast(_device_data->_d_py.data()),
             thrust::raw_pointer_cast(_device_data->_d_pz.data()),
             thrust::raw_pointer_cast(_device_data->_d_fx.data()),
             thrust::raw_pointer_cast(_device_data->_d_fy.data()),
             thrust::raw_pointer_cast(_device_data->_d_fz.data()),
             thrust::raw_pointer_cast(_device_data->_d_evdwl.data()));

    // // 计算总势能
    //rbmd::Real* d_total_evdwl;
    //hipMalloc(&d_total_evdwl, sizeof(rbmd::Real));
    //op::DeviceReduceSum(thrust::raw_pointer_cast(_device_data->_d_evdwl.data()), d_total_evdwl, _num_atoms);

    //// 从设备复制结果到主机
    //rbmd::Real total_evdwl;
    //hipMemcpy(&total_evdwl, d_total_evdwl, sizeof(rbmd::Real), hipMemcpyDeviceToHost);

    //std::cout << "Total evdwl: " << total_evdwl << std::endl;

    //hipFree(d_total_evdwl);

             
}

