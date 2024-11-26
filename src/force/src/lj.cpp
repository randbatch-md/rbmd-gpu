#include "lj.h"

#include <thrust/device_ptr.h>

#include "../../common/device_types.h"
#include "../../common/rbmd_define.h"
#include "../../common/types.h"
#include "../data_manager/include/model/md_data.h"
#include "lj_op/lj_op.h"
#include "neighbor_list/include/neighbor_list_builder/full_neighbor_list_builder.h"
#include "neighbor_list/include/neighbor_list_builder/half_neighbor_list_builder.h"
#include "neighbor_list/include/neighbor_list_builder/rbl_full_neighbor_list_builder.h"
// #include <hipcub/hipcub.hpp>
// #include <hipcub/backend/rocprim/block/block_reduce.hpp>
extern int test_current_step;
rbmd::Real test_ave_pe_rbl;
rbmd::Real test_ave_pe_init;
LJ::LJ() {
  _rbl_neighbor_list_builder = std::make_shared<RblFullNeighborListBuilder>();
  _neighbor_list_builder = std::make_shared<FullNeighborListBuilder>();
  std::remove("thermo_local.txt");
}

LJ::~LJ()
{
}

void LJ::Init() {
  _cut_off = DataManager::getInstance().getConfigData()->Get
 <rbmd::Real>("cut_off", "hyper_parameters", "neighbor");

  _neighbor_type =
    DataManager::getInstance().getConfigData()->Get<std::string>(
        "type", "hyper_parameters", "neighbor");
}

void LJ::Execute()
{
  if (_neighbor_type == "RBL")  // RBL
  {
    ComputeLJRBL();
  }
  else  //
  {
    ComputeLJVerlet();
  }

  //
  EvaluatePotentialenergy();
}

void LJ::ComputeLJRBL()
{
    // rbl_neighbor_list_build
    auto start = std::chrono::high_resolution_clock::now();
    _rbl_list = _rbl_neighbor_list_builder->Build();

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<rbmd::Real> duration = end - start;
    std::cout << "构建RBL邻居列表耗时" << duration.count() << "秒" << std::endl;

    // compute force
    const auto r_core =
        DataManager::getInstance().getConfigData()->Get<rbmd::Real>(
            "r_core", "hyper_parameters", "neighbor");

    const auto neighbor_sample_num =
        DataManager::getInstance().getConfigData()->Get<rbmd::Id>(
            "neighbor_sample_num", "hyper_parameters", "neighbor");

    auto num_atoms = *(_structure_info_data->_num_atoms);
    op::LJRBLForceOp<device::DEVICE_GPU>()(
        *_box, r_core, _cut_off,
        num_atoms,neighbor_sample_num,_rbl_list->_selection_frequency,
        thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
        thrust::raw_pointer_cast(_device_data->_d_sigma.data()),
        thrust::raw_pointer_cast(_device_data->_d_eps.data()),
        thrust::raw_pointer_cast(_rbl_list->_start_idx.data()),
        thrust::raw_pointer_cast(_rbl_list->_end_idx.data()),
        thrust::raw_pointer_cast(_rbl_list->_d_neighbors.data()),
        thrust::raw_pointer_cast(_rbl_list->_d_random_neighbor.data()),
        thrust::raw_pointer_cast(_rbl_list->_d_random_neighbor_num.data()),
        thrust::raw_pointer_cast(_device_data->_d_px.data()),
        thrust::raw_pointer_cast(_device_data->_d_py.data()),
        thrust::raw_pointer_cast(_device_data->_d_pz.data()),
        thrust::raw_pointer_cast(_device_data->_d_fx.data()),
        thrust::raw_pointer_cast(_device_data->_d_fy.data()),
        thrust::raw_pointer_cast(_device_data->_d_fz.data()));

    _corr_value_x =
        thrust::reduce(_device_data->_d_fx.begin(), _device_data->_d_fx.end(),
                       0.0f, thrust::plus<rbmd::Real>()) /num_atoms;
    _corr_value_y =
        thrust::reduce(_device_data->_d_fy.begin(), _device_data->_d_fy.end(),
                       0.0f, thrust::plus<rbmd::Real>()) /num_atoms;
    _corr_value_z =
        thrust::reduce(_device_data->_d_fz.begin(), _device_data->_d_fz.end(),
                       0.0f, thrust::plus<rbmd::Real>()) /num_atoms;

    // fix RBL:   rbl_force = force - corr_value
    op::FixRBLForceOp<device::DEVICE_GPU>()(
                        num_atoms, _corr_value_x, _corr_value_y, _corr_value_z,
                        thrust::raw_pointer_cast(_device_data->_d_fx.data()),
                        thrust::raw_pointer_cast(_device_data->_d_fy.data()),
                        thrust::raw_pointer_cast(_device_data->_d_fz.data()));

    //energy
    ComputeLJEnergy();
}

void LJ::ComputeLJVerlet()
{
  // neighbor_list_build
  auto start = std::chrono::high_resolution_clock::now();
  _list = _neighbor_list_builder->Build();

  auto end = std::chrono::high_resolution_clock::now();

  std::chrono::duration<rbmd::Real> duration = end - start;
  std::cout << "构建verlet-list耗时" << duration.count() << "秒" << std::endl;

  //
  thrust::device_vector<rbmd::Real> d_total_evdwl(1, 0.0);
  auto num_atoms = *(_structure_info_data->_num_atoms);
  // compute LJ
  op::LJForceOp<device::DEVICE_GPU>()(
              *_box, _cut_off,num_atoms,
              thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
              thrust::raw_pointer_cast(_device_data->_d_sigma.data()),
              thrust::raw_pointer_cast(_device_data->_d_eps.data()),
              thrust::raw_pointer_cast(_list->_start_idx.data()),
              thrust::raw_pointer_cast(_list->_end_idx.data()),
              thrust::raw_pointer_cast(_list->_d_neighbors.data()),
              thrust::raw_pointer_cast(_device_data->_d_px.data()),
              thrust::raw_pointer_cast(_device_data->_d_py.data()),
              thrust::raw_pointer_cast(_device_data->_d_pz.data()),
              thrust::raw_pointer_cast(_device_data->_d_fx.data()),
              thrust::raw_pointer_cast(_device_data->_d_fy.data()),
              thrust::raw_pointer_cast(_device_data->_d_fz.data()),
              thrust::raw_pointer_cast(_device_data->_d_flat_virial.data()),
              thrust::raw_pointer_cast(d_total_evdwl.data()));

  // 从设备端拷贝数据到主机端
  thrust::host_vector<rbmd::Real> h_total_evdwl(d_total_evdwl);
  _ave_evdwl = h_total_evdwl[0] / num_atoms;

  std::cout << "test_current_step:" << test_current_step << " "
            << "average_vdwl_energy:" << _ave_evdwl << std::endl;
  std::cout << "out of force execute" << std::endl;

  //sum virial on host
  std::vector<rbmd::Real> h_total_virial(num_atoms * 6);
  thrust::copy(_device_data->_d_flat_virial.begin(),
    _device_data->_d_flat_virial.end(), h_total_virial.begin());

  std::vector<rbmd::Real> virial(6);
  virial =  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  for(int atom = 0; atom < num_atoms; ++atom){
    for(int i = 0; i < 6; ++i){
      virial[i] += h_total_virial[atom * 6 + i];
    }
  }

  thrust::copy(virial.begin(),
  virial.end(), _device_data->_d_virial_lj.begin());

}

void LJ::ComputeLJEnergy()
{
  // energy
  _list = _neighbor_list_builder->Build();

  thrust::device_vector<rbmd::Real> d_total_evdwl(1, 0.0);
  auto num_atoms = *(_structure_info_data->_num_atoms);
  op::LJEnergyOp<device::DEVICE_GPU>()(
                *_box, _cut_off, num_atoms,
               thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
               thrust::raw_pointer_cast(_device_data->_d_sigma.data()),
               thrust::raw_pointer_cast(_device_data->_d_eps.data()),
               thrust::raw_pointer_cast(_list->_start_idx.data()),
               thrust::raw_pointer_cast(_list->_end_idx.data()),
               thrust::raw_pointer_cast(_list->_d_neighbors.data()),
               thrust::raw_pointer_cast(_device_data->_d_px.data()),
               thrust::raw_pointer_cast(_device_data->_d_py.data()),
               thrust::raw_pointer_cast(_device_data->_d_pz.data()),
               thrust::raw_pointer_cast(_device_data->_d_flat_virial.data()),
               thrust::raw_pointer_cast(d_total_evdwl.data()));

  // 从设备端拷贝数据到主机端
  thrust::host_vector<rbmd::Real> h_total_evdwl(d_total_evdwl);
  _ave_evdwl = h_total_evdwl[0] / num_atoms;

  std::cout << "test_current_step:" << test_current_step << " "
            << "average_vdwl_energy:" << _ave_evdwl << std::endl;

  //sum virial on host
  std::vector<rbmd::Real> h_total_virial(num_atoms * 6);
  thrust::copy(_device_data->_d_flat_virial.begin(),
    _device_data->_d_flat_virial.end(), h_total_virial.begin());

  std::vector<rbmd::Real> virial(6);
  virial =  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  for(int atom = 0; atom < num_atoms; ++atom){
    for(int i = 0; i < 6; ++i){
      virial[i] += h_total_virial[atom * 6 + i];
    }
  }


  thrust::copy(virial.begin(),
    virial.end(), _device_data->_d_virial_lj.begin());

}

void LJ::EvaluatePotentialenergy()
{
  _ave_pe_rbl = _ave_evdwl_rbl;
  test_ave_pe_rbl = _ave_pe_rbl;


  if(1 == test_current_step)
  {
    _ave_pe_init = _ave_evdwl;
    test_ave_pe_init = _ave_pe_init;
  }
  _ave_pe = _ave_evdwl;

  //out
  std::ofstream outfile("thermo_local.txt", std::ios::app);
  if (outfile.tellp() == 0) {
    outfile << "step _ave_pe_rbl  _ave_pe" << std::endl;
  }
  outfile << test_current_step << " " << _ave_pe_rbl  << " "<< _ave_pe << std::endl;
  outfile.close();
}


