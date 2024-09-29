#include "ljforce.h"

#include <thrust/device_ptr.h>

#include "../../common/device_types.h"
#include "../../common/rbmd_define.h"
#include "../../common/types.h"
#include "ljforce_op/ljforce_op.h"
#include "neighbor_list/include/neighbor_list_builder/half_neighbor_list_builder.h"
#include "neighbor_list/include/neighbor_list_builder/full_neighbor_list_builder.h"
#include "neighbor_list/include/neighbor_list_builder/rbl_full_neighbor_list_builder.h"
// #include <hipcub/hipcub.hpp>
// #include <hipcub/backend/rocprim/block/block_reduce.hpp>

LJForce::LJForce()
{
  _rbl_neighbor_list_builder = std::make_shared<RblFullNeighborListBuilder>();
  _neighbor_list_builder = std::make_shared<FullNeighborListBuilder>();

  CHECK_RUNTIME(MALLOC(&_d_total_evdwl, sizeof(rbmd::Real)));

}

LJForce::~LJForce()
{
    CHECK_RUNTIME(FREE(_d_total_evdwl));
}

void LJForce::Init() 
{ 
    _num_atoms = *(_structure_info_data->_num_atoms);
    _corr_value_x = 0;
    _corr_value_y = 0;
    _corr_value_z = 0;
}

void LJForce::Execute() 
{
  extern int test_current_step;
  const auto cut_off = DataManager::getInstance().getConfigData()->Get
      <rbmd::Real>("cut_off", "hyper_parameters", "neighbor");

  auto neighbor_type = DataManager::getInstance().getConfigData()->Get
      <std::string>("type", "hyper_parameters", "neighbor");
  if (neighbor_type == "RBL")   //RBL
  {
      //rbl_neighbor_list_build
      auto start = std::chrono::high_resolution_clock::now();
      rbl_list = _rbl_neighbor_list_builder->Build();

      auto end = std::chrono::high_resolution_clock::now();

      std::chrono::duration<rbmd::Real> duration = end - start;
      std::cout << "构建RBL邻居列表耗时" << duration.count() << "秒" << std::endl;


      //compute force
      const auto r_core  = DataManager::getInstance().getConfigData()->Get
          <rbmd::Real>( "r_core", "hyper_parameters", "neighbor");

      const auto neighbor_sample_num = DataManager::getInstance().getConfigData()->Get<rbmd::Id>("neighbor_sample_num", "hyper_parameters", "neighbor");

      op::LJRBLForceOp<device::DEVICE_GPU> lj_rbl_force_op;
      lj_rbl_force_op(_device_data->_d_box, r_core, cut_off, _num_atoms, neighbor_sample_num, rbl_list->_selection_frequency,
          thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
          thrust::raw_pointer_cast(_device_data->_d_molecular_id.data()),
          thrust::raw_pointer_cast(_device_data->_d_sigma.data()),
          thrust::raw_pointer_cast(_device_data->_d_eps.data()),
          thrust::raw_pointer_cast(rbl_list->_start_idx.data()),
          thrust::raw_pointer_cast(rbl_list->_end_idx.data()),
          thrust::raw_pointer_cast(rbl_list->_d_neighbors.data()),
          thrust::raw_pointer_cast(rbl_list->_d_random_neighbor.data()),
          thrust::raw_pointer_cast(rbl_list->_d_random_neighbor_num.data()),
          thrust::raw_pointer_cast(_device_data->_d_px.data()),
          thrust::raw_pointer_cast(_device_data->_d_py.data()),
          thrust::raw_pointer_cast(_device_data->_d_pz.data()),
          thrust::raw_pointer_cast(_device_data->_d_fx.data()),
          thrust::raw_pointer_cast(_device_data->_d_fy.data()),
          thrust::raw_pointer_cast(_device_data->_d_fz.data()));

      _corr_value_x = thrust::reduce(_device_data->_d_fx.begin(), _device_data->_d_fx.end(), 0.0f, thrust::plus<rbmd::Real>()) / _num_atoms;
      _corr_value_y = thrust::reduce(_device_data->_d_fy.begin(), _device_data->_d_fy.end(), 0.0f, thrust::plus<rbmd::Real>()) / _num_atoms;
      _corr_value_z = thrust::reduce(_device_data->_d_fz.begin(), _device_data->_d_fz.end(), 0.0f, thrust::plus<rbmd::Real>()) / _num_atoms;


      //fix RBL:   rbl_force = f - corr_value
      op::FixLJRBLForceOp<device::DEVICE_GPU> fix_lj_rbl_force_op;
      fix_lj_rbl_force_op(_num_atoms, _corr_value_x, _corr_value_y, _corr_value_z,
          thrust::raw_pointer_cast(_device_data->_d_fx.data()),
          thrust::raw_pointer_cast(_device_data->_d_fy.data()),
          thrust::raw_pointer_cast(_device_data->_d_fz.data()));


      //energy
      list = _neighbor_list_builder->Build();

      rbmd::Real h_total_evdwl = 0.0;

      CHECK_RUNTIME(MEMSET(_d_total_evdwl, 0, sizeof(rbmd::Real)));

      op::LJEnergyOp<device::DEVICE_GPU> lj_energy_op;
      lj_energy_op(_device_data->_d_box, cut_off, _num_atoms,
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
          _d_total_evdwl);

      CHECK_RUNTIME(MEMCPY(&h_total_evdwl, _d_total_evdwl, sizeof(rbmd::Real), D2H));

      // 打印累加后的总能量
      rbmd::Real  ave_evdwl = h_total_evdwl / _num_atoms;
      std::cout << "test_current_step:" << test_current_step << " " << "average_vdwl_energy:" << ave_evdwl << std::endl;

      ////out
      std::ofstream outfile("ave_evdwl.txt", std::ios::app);
      outfile << test_current_step << " " << ave_evdwl << std::endl;
      outfile.close();
  }

  else  //
  {
      //neighbor_list_build
      auto start = std::chrono::high_resolution_clock::now();
      list = _neighbor_list_builder->Build();

      auto end = std::chrono::high_resolution_clock::now();

      std::chrono::duration<rbmd::Real> duration = end - start;
      std::cout << "构建verlet-list耗时" << duration.count() << "秒" << std::endl;

      //

      rbmd::Real h_total_evdwl = 0.0;

      CHECK_RUNTIME(MEMSET(_d_total_evdwl, 0, sizeof(rbmd::Real)));

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
                  thrust::raw_pointer_cast(_device_data->_d_virial_xx.data()),
                  thrust::raw_pointer_cast(_device_data->_d_virial_yy.data()),
                  thrust::raw_pointer_cast(_device_data->_d_virial_zz.data()),
                  thrust::raw_pointer_cast(_device_data->_d_virial_xy.data()),
                  thrust::raw_pointer_cast(_device_data->_d_virial_xz.data()),
                  thrust::raw_pointer_cast(_device_data->_d_virial_yz.data()),                
                  _d_total_evdwl);

      virial[0] = thrust::reduce(_device_data->_d_virial_xx.begin(), _device_data->_d_virial_xx.end(), 0.0f, thrust::plus<rbmd::Real>());
      virial[1] = thrust::reduce(_device_data->_d_virial_yy.begin(), _device_data->_d_virial_yy.end(), 0.0f, thrust::plus<rbmd::Real>());
      virial[2] = thrust::reduce(_device_data->_d_virial_zz.begin(), _device_data->_d_virial_zz.end(), 0.0f, thrust::plus<rbmd::Real>());
      virial[3] = thrust::reduce(_device_data->_d_virial_xy.begin(), _device_data->_d_virial_xy.end(), 0.0f, thrust::plus<rbmd::Real>());
      virial[4] = thrust::reduce(_device_data->_d_virial_xz.begin(), _device_data->_d_virial_xz.end(), 0.0f, thrust::plus<rbmd::Real>());
      virial[5] = thrust::reduce(_device_data->_d_virial_yz.begin(), _device_data->_d_virial_yz.end(), 0.0f, thrust::plus<rbmd::Real>());

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
}

void LJForce::ComputeChargeStructureFactorEwald(
    Box* box,
    rbmd::Id num_atoms,
    rbmd::Id Kmax, 
    rbmd::Real* value_Re_array,
    rbmd::Real* value_Im_array)
{
    thrust::device_vector<rbmd::Real> density_real(num_atoms);
    thrust::device_vector<rbmd::Real> density_imag(num_atoms);
  
    thrust::fill(density_real.begin(), density_real.end(), 0.0f);
    thrust::fill(density_imag.begin(), density_imag.end(), 0.0f);
    
    rbmd::Id total_elements = (2 * Kmax + 1) * (2 * Kmax + 1) * (2 * Kmax + 1) - 1;

    rbmd::Id index = 0;
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

                    value_Re_array[index] = value_Re;
                    value_Im_array[index] = value_Im;
                    index++; // 增加索引
                }
            }
        }
    }
}

void LJForce::ComputeEwladForce() 
{
    rbmd::Id  Kmax = 2;
    rbmd::Real  alpha = 0.01;
    rbmd::Real* value_Re_array;
    rbmd::Real* value_Im_array;
    rbmd::Id total_K_elements = (2 * Kmax + 1) * (2 * Kmax + 1) * (2 * Kmax + 1) - 1;

    CHECK_RUNTIME(MALLOC(&value_Re_array, total_K_elements * sizeof(rbmd::Real)));
    CHECK_RUNTIME(MALLOC(&value_Im_array, total_K_elements * sizeof(rbmd::Real)));


    ComputeChargeStructureFactorEwald(_device_data->_d_box, _num_atoms, Kmax, value_Re_array, value_Im_array);

    op::ComputeEwaldForceOp<device::DEVICE_GPU> ewlad_force_op;
    ewlad_force_op(_device_data->_d_box,_num_atoms, Kmax, alpha,
        value_Re_array,
        value_Im_array,
        thrust::raw_pointer_cast(_device_data->_d_charge.data()),
        thrust::raw_pointer_cast(_device_data->_d_px.data()),
        thrust::raw_pointer_cast(_device_data->_d_py.data()),
        thrust::raw_pointer_cast(_device_data->_d_pz.data()),
        thrust::raw_pointer_cast(_device_data->_d_fx.data()),
        thrust::raw_pointer_cast(_device_data->_d_fy.data()),
        thrust::raw_pointer_cast(_device_data->_d_fz.data()));

    CHECK_RUNTIME(FREE(value_Re_array));
    CHECK_RUNTIME(FREE(value_Im_array));

}

