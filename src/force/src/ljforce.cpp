#include "ljforce.h"
#include <thrust/device_ptr.h>
#include "../../common/device_types.h"
#include "../../common/types.h"
#include "ljforce_op/ljforce_op.h"
#include <thrust/reduce.h>
// #include <hipcub/hipcub.hpp>
// #include <hipcub/backend/rocprim/block/block_reduce.hpp>

LJForce::LJForce()
{
  full_list_builder = std::make_shared<FullNeighborListBuilder>();
};

void LJForce::Init() { _num_atoms = *(_structure_info_data->_num_atoms);}

void LJForce::Execute() 
{

  rbmd::Real cut_off = 5.0;
  //
  extern int test_current_step;
  list = full_list_builder->Build();
  //list->print("./neighbor_list_step_" + std::to_string(test_current_step)+"_.csv");

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
  std::cout << "test_current_step:" << test_current_step <<  " " << "average_vdwl_energy:" << ave_evdwl << std::endl;

  // 释放设备端分配的内存
  CHECK_RUNTIME(FREE(d_total_evdwl));

  std::cout << "out of force execute" << std::endl;


  //out
  std::ofstream outfile("ave_evdwl.txt", std::ios::app);
  outfile << test_current_step << " " << ave_evdwl << std::endl;
  outfile.close();   
}

void LJForce::ComputeChargeStructureFactorEwald(
    Box* box,
    rbmd::Id num_atoms,
    rbmd::Id Kmax, 
    std::vector<float2> rhok)
{

    thrust::device_vector<rbmd::Real> density_real(num_atoms, 0.0f);
    thrust::device_vector<rbmd::Real> density_imag(num_atoms, 0.0f);

    for (rbmd::Id i = -Kmax; i <= Kmax; i++)
    {
        for (rbmd::Id j = -Kmax; j <= Kmax; j++)
        {
            for (rbmd::Id k = -Kmax; k <= Kmax; k++)
            {
                if (!(i == 0 && j == 0 && k == 0))
                {
                    float3 K = make_float3(2.0 * M_PI * i / box->_length[0],
                                           2.0 * M_PI * j / box->_length[1],
                                           2.0 * M_PI * k / box->_length[2]);

                    op::ComputeChargeStructureFactorComponentOp<device::DEVICE_GPU> charge_structure_factor_op;

                    charge_structure_factor_op(num_atoms, K,
                        thrust::raw_pointer_cast(_device_data->_d_px.data()),
                        thrust::raw_pointer_cast(_device_data->_d_py.data()),
                        thrust::raw_pointer_cast(_device_data->_d_pz.data()),
                        thrust::raw_pointer_cast(_device_data->_d_charge.data()),
                        thrust::raw_pointer_cast(density_real.data()),
                        thrust::raw_pointer_cast(density_imag.data()));

                    rbmd::Real value_Re = thrust::reduce(density_real.begin(), density_real.end(), 0.0f, thrust::plus<rbmd::Real>());
                    rbmd::Real value_Im = thrust::reduce(density_imag.begin(), density_imag.end(), 0.0f, thrust::plus<rbmd::Real>());

                    rhok.push_back(make_float2(value_Re, value_Im));
                }
            }
        }
    }
}

