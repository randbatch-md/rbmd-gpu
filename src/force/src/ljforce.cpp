#include "ljforce.h"

#include <thrust/device_ptr.h>

#include "../../common/device_types.h"
#include "../../common/rbmd_define.h"
#include "../../common/types.h"
#include "ljforce_op/ljforce_op.h"
#include "neighbor_list/include/neighbor_list_builder/full_neighbor_list_builder.h"
#include "neighbor_list/include/neighbor_list_builder/half_neighbor_list_builder.h"
#include "neighbor_list/include/neighbor_list_builder/rbl_full_neighbor_list_builder.h"
// #include <hipcub/hipcub.hpp>
// #include <hipcub/backend/rocprim/block/block_reduce.hpp>
extern int test_current_step;

LJForce::LJForce() {
  _rbl_neighbor_list_builder = std::make_shared<RblFullNeighborListBuilder>();
  _neighbor_list_builder = std::make_shared<FullNeighborListBuilder>();

  CHECK_RUNTIME(MALLOC(&_d_total_evdwl, sizeof(rbmd::Real)));
}

LJForce::~LJForce()
{
  CHECK_RUNTIME(FREE(_d_total_evdwl));
}

void LJForce::Init() {
  _num_atoms = *(_structure_info_data->_num_atoms);
  _corr_value_x = 0;
  _corr_value_y = 0;
  _corr_value_z = 0;

  _cut_off = DataManager::getInstance().getConfigData()->Get
 <rbmd::Real>("cut_off", "hyper_parameters", "neighbor");

  _neighbor_type =
    DataManager::getInstance().getConfigData()->Get<std::string>(
        "type", "hyper_parameters", "neighbor");
}

void LJForce::Execute() {


  if (_neighbor_type == "RBL")  // RBL
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

    op::LJRBLForceOp<device::DEVICE_GPU> lj_rbl_force_op;
    lj_rbl_force_op(
        _device_data->_d_box, r_core, _cut_off, _num_atoms, neighbor_sample_num,
        _rbl_list->_selection_frequency,
        thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
        thrust::raw_pointer_cast(_device_data->_d_molecular_id.data()),
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
                       0.0f, thrust::plus<rbmd::Real>()) /_num_atoms;
    _corr_value_y =
        thrust::reduce(_device_data->_d_fy.begin(), _device_data->_d_fy.end(),
                       0.0f, thrust::plus<rbmd::Real>()) /_num_atoms;
    _corr_value_z =
        thrust::reduce(_device_data->_d_fz.begin(), _device_data->_d_fz.end(),
                       0.0f, thrust::plus<rbmd::Real>()) /_num_atoms;

    // fix RBL:   rbl_force = f - corr_value
    op::FixRBLForceOp<device::DEVICE_GPU> fix_rbl_force_op;
    fix_rbl_force_op(_num_atoms, _corr_value_x, _corr_value_y, _corr_value_z,
                        thrust::raw_pointer_cast(_device_data->_d_fx.data()),
                        thrust::raw_pointer_cast(_device_data->_d_fy.data()),
                        thrust::raw_pointer_cast(_device_data->_d_fz.data()));

    //energy
    //ComputeLJEnergy();
  }

  else  //
  {
    // neighbor_list_build
    auto start = std::chrono::high_resolution_clock::now();
    _list = _neighbor_list_builder->Build();

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<rbmd::Real> duration = end - start;
    std::cout << "构建verlet-list耗时" << duration.count() << "秒" << std::endl;

    //

    rbmd::Real h_total_evdwl = 0.0;

    CHECK_RUNTIME(MEMSET(_d_total_evdwl, 0, sizeof(rbmd::Real)));

    // compute LJForce
    op::LJForceOp<device::DEVICE_GPU> lj_force_op;
    lj_force_op(_device_data->_d_box, _cut_off, _num_atoms,
                thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
                thrust::raw_pointer_cast(_device_data->_d_molecular_id.data()),
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
                _d_total_evdwl);

    // CHECK_RUNTIME(
    //     MEMCPY(&h_total_evdwl, _d_total_evdwl, sizeof(rbmd::Real), D2H));
    //
    // // 打印累加后的总能量
    // rbmd::Real ave_evdwl = h_total_evdwl / _num_atoms;
    // std::cout << "test_current_step:" << test_current_step << " "
    //           << "average_vdwl_energy:" << ave_evdwl << std::endl;
    //
    // std::cout << "out of force execute" << std::endl;
    //
    // // out
    // std::ofstream outfile("ave_evdwl.txt", std::ios::app);
    // outfile << test_current_step << " " << ave_evdwl << std::endl;
    // outfile.close();
  }
}

void LJForce::ComputeLJEnergy()
{
  // energy
  _list = _neighbor_list_builder->Build();

  rbmd::Real h_total_evdwl = 0.0;

  CHECK_RUNTIME(MEMSET(_d_total_evdwl, 0, sizeof(rbmd::Real)));

  op::LJEnergyOp<device::DEVICE_GPU> lj_energy_op;
  lj_energy_op(_device_data->_d_box, _cut_off, _num_atoms,
               thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
               thrust::raw_pointer_cast(_device_data->_d_molecular_id.data()),
               thrust::raw_pointer_cast(_device_data->_d_sigma.data()),
               thrust::raw_pointer_cast(_device_data->_d_eps.data()),
               thrust::raw_pointer_cast(_list->_start_idx.data()),
               thrust::raw_pointer_cast(_list->_end_idx.data()),
               thrust::raw_pointer_cast(_list->_d_neighbors.data()),
               thrust::raw_pointer_cast(_device_data->_d_px.data()),
               thrust::raw_pointer_cast(_device_data->_d_py.data()),
               thrust::raw_pointer_cast(_device_data->_d_pz.data()),
               _d_total_evdwl);

  CHECK_RUNTIME(
      MEMCPY(&h_total_evdwl, _d_total_evdwl, sizeof(rbmd::Real), D2H));

  // 打印累加后的总能量
  rbmd::Real ave_evdwl = h_total_evdwl / _num_atoms;
  std::cout << "test_current_step:" << test_current_step << " "
            << "average_vdwl_energy:" << ave_evdwl << std::endl;

  ////out
  std::ofstream outfile("ave_evdwl.txt", std::ios::app);
  outfile << test_current_step << " " << ave_evdwl << std::endl;
  outfile.close();
}


