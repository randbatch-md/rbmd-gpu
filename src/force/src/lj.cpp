#include "lj.h"

#include <output/include/Logger.hpp>

#include "../../common/device_types.h"
#include "../../common/rbmd_define.h"
#include "../../common/types.h"
#include "lj_op/lj_op.h"
#include "neighbor_list/include/neighbor_list_builder/full_neighbor_list_builder.h"
#include "neighbor_list/include/neighbor_list_builder/rbl_full_neighbor_list_builder.h"
#include "common/timing_statistics.hpp"
#include "common/thermo_stats.hpp"
// #include <hipcub/hipcub.hpp>
// #include <hipcub/backend/rocprim/block/block_reduce.hpp>
extern int test_current_step;
rbmd::Real test_e_pe_rbl;
rbmd::Real test_e_pe_init;
LJ::LJ() {
  _rbl_neighbor_list_builder = std::make_shared<RblFullNeighborListBuilder>();
  _neighbor_list_builder = std::make_shared<FullNeighborListBuilder>();
  std::remove("thermo.txt");
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
    TimingStatistics::Instance().record("Neighbor-List",duration.count());

    // compute force
    auto start_rbl_force = std::chrono::high_resolution_clock::now();
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

  auto end_rbl_force = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration_rbl_force = end_rbl_force - start_rbl_force;
  TimingStatistics::Instance().record("Short-Range",duration_rbl_force.count());

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
  TimingStatistics::Instance().record("Neighbor-List",duration.count());
  //
  auto start_verlet_force = std::chrono::high_resolution_clock::now();
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

  auto end_verlet_force = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration_verlet_force = end_verlet_force - start_verlet_force;
  TimingStatistics::Instance().record("Short-Range",duration_verlet_force.count());

  // D2H
  thrust::host_vector<rbmd::Real> h_total_evdwl(d_total_evdwl);
  _e_vdwl = h_total_evdwl[0] / num_atoms;

  ThermoStats::Instance().AddThermoData("vdwl",_e_vdwl);

  //sum virial_lj on host
  ReduceVirial(num_atoms,_device_data->_d_flat_virial,
_device_data->_d_virial_lj);
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

  // D2H
  thrust::host_vector<rbmd::Real> h_total_evdwl(d_total_evdwl);
  _e_vdwl = h_total_evdwl[0] / num_atoms;

  ThermoStats::Instance().AddThermoData("vdwl",_e_vdwl);


  //sum virial_lj on host
  ReduceVirial(num_atoms,_device_data->_d_flat_virial,
_device_data->_d_virial_lj);
}

void LJ::EvaluatePotentialenergy()
{
  _e_pe_rbl = _e_vdwl_rbl;
  test_e_pe_rbl = _e_pe_rbl;


  if(1 == test_current_step)
  {
    _e_pe_init = _e_vdwl;
    test_e_pe_init = _e_pe_init;
  }
  _e_pe = _e_vdwl;

  //out
  std::ofstream outfile("thermo.txt", std::ios::app);
  if (outfile.tellp() == 0) {
    outfile << "step e_vdwl  e_pe" << std::endl;
  }
  outfile << test_current_step << " " << _e_vdwl  << " "<< _e_pe << std::endl;
  outfile.close();
}


