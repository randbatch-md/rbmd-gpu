#include "lj_cut_coul_kspace.h"

#include "../../common/device_types.h"
#include "../../common/rbmd_define.h"
#include "../../common/types.h"
#include "../../common/unit_factor.h"
#include "force_op/force_op.h"
#include "lj_op/lj_op.h"
#include "lj_cut_coul_kspace_op/lj_cut_coul_kspace_op.h"
#include "../common/RBEPSample.h"
#include "neighbor_list/include/linked_cell/linked_cell_locator.h"
#include "neighbor_list/include/neighbor_list_builder/full_neighbor_list_builder.h"
#include "neighbor_list/include/neighbor_list_builder/rbl_full_neighbor_list_builder.h"
// #include <hipcub/hipcub.hpp>
// #include <hipcub/backend/rocprim/block/block_reduce.hpp>
#include "common/thermo_stats.hpp"
#include "common/timing_statistics.hpp"


extern int test_current_step;
extern std::map<std::string, UNIT> unit_factor_map;

LJCutCoulKspace::LJCutCoulKspace()
{
  _rbl_neighbor_list_builder = std::make_shared<RblFullNeighborListBuilder>();
  _neighbor_list_builder = std::make_shared<FullNeighborListBuilder>();

  // 创建 K-Space 计算器实例，并将自身所需的数据和指针传进去
  _kspace_calculator = std::make_unique<KSpaceCalculator>();

  auto unit = DataManager::getInstance().getConfigData()->Get
<std::string>("unit", "init_configuration", "read_data");
  UNIT unit_factor = unit_factor_map[unit];
  switch (unit_factor) {
    case UNIT::LJ:
      _qqr2e = UnitFactor<UNIT::LJ>::_qqr2e;
      break;
    case UNIT::METAL:
      _qqr2e = UnitFactor<UNIT::METAL>::_qqr2e;
      break;
    case UNIT::REAL:
      _qqr2e = UnitFactor<UNIT::REAL>::_qqr2e;
      break;

    default:
      break;
  }

  std::remove("thermo.txt");
}

LJCutCoulKspace::~LJCutCoulKspace(){}

void LJCutCoulKspace::Init()
{
  const auto& config = DataManager::getInstance().getConfigData();

  //neighbor
  _cut_off = config->Get<rbmd::Real>("cut_off", "hyper_parameters", "neighbor");
  _neighbor_type = config->Get<std::string>("type", "hyper_parameters", "neighbor");

  if("RBL" == _neighbor_type) {
    bool energy_rbl_flag = config->PathExists({"hyper_parameters", "neighbor" ,"energy_rbl_flag"});
    if (energy_rbl_flag) {
      _energy_rbl_flag = config->Get<std::string>("energy_rbl_flag", "hyper_parameters", "neighbor");
    }
    else {
      Logger::Instance().error( "\033[31m When using RBL for the neighbor type, "
                   "the key 'energy_rbl_flag' must be defined.\033[0m");
      exit(EXIT_FAILURE); //
    }
  }

  // 2. 初始化 K-Space 组件
  _alpha = config->Get<rbmd::Real>("alpha", "hyper_parameters", "coulomb");
  _kspace_calculator->Init();
}

void LJCutCoulKspace::Execute()
{
  ComputeLJCutCoulForce();
  // 2. 委托 K-Space 组件计算长程力
  _kspace_calculator->Execute();
  SumForces();

  EvaluatePotentialenergy();
}

void LJCutCoulKspace::ComputeLJCutCoulForce()
{
  //
  if ("RBL" ==_neighbor_type)
  {
    ComputeLJRBL();
  }
  else
  {
    ComputeLJVerlet();
  }

  //add thermo
  ThermoStats::Instance().AddThermoData("vdwl",_e_vdwl);
  ThermoStats::Instance().AddThermoData("coul",_e_coul);
}

void LJCutCoulKspace::ComputeLJRBL()
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
    op::LJCutCoulRBLForceOp<device::DEVICE_GPU>()(
        *_box,_device_data->_d_erf_table,
        r_core, _cut_off,num_atoms,neighbor_sample_num,
        _rbl_list->_selection_frequency,_alpha,_qqr2e,
        thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
        thrust::raw_pointer_cast(_device_data->_d_sigma.data()),
        thrust::raw_pointer_cast(_device_data->_d_eps.data()),
        thrust::raw_pointer_cast(_rbl_list->_start_idx.data()),
        thrust::raw_pointer_cast(_rbl_list->_end_idx.data()),
        thrust::raw_pointer_cast(_rbl_list->_d_neighbors.data()),
        thrust::raw_pointer_cast(_rbl_list->_d_random_neighbor.data()),
        thrust::raw_pointer_cast(_rbl_list->_d_random_neighbor_num.data()),
        thrust::raw_pointer_cast(_device_data->_d_charge.data()),
        thrust::raw_pointer_cast(_device_data->_d_px.data()),
        thrust::raw_pointer_cast(_device_data->_d_py.data()),
        thrust::raw_pointer_cast(_device_data->_d_pz.data()),
        thrust::raw_pointer_cast(_device_data->_d_force_ljcoul_x.data()),
        thrust::raw_pointer_cast(_device_data->_d_force_ljcoul_y.data()),
        thrust::raw_pointer_cast(_device_data->_d_force_ljcoul_z.data()));

    _corr_value_x =
        thrust::reduce(_device_data->_d_force_ljcoul_x.begin(),
          _device_data->_d_force_ljcoul_x.end(),
                         0.0f, thrust::plus<rbmd::Real>()) /num_atoms;
    _corr_value_y =
        thrust::reduce(_device_data->_d_force_ljcoul_y.begin(),
          _device_data->_d_force_ljcoul_y.end(),
                       0.0f, thrust::plus<rbmd::Real>()) /num_atoms;
    _corr_value_z =
        thrust::reduce(_device_data->_d_force_ljcoul_z.begin(),
          _device_data->_d_force_ljcoul_z.end(),
                       0.0f, thrust::plus<rbmd::Real>()) /num_atoms;

    // fix RBL:   rbl_force = force - corr_value
    op::FixRBLForceOp<device::DEVICE_GPU>()(
                         num_atoms, _corr_value_x, _corr_value_y, _corr_value_z,
                        thrust::raw_pointer_cast(_device_data->_d_force_ljcoul_x.data()),
                        thrust::raw_pointer_cast(_device_data->_d_force_ljcoul_y.data()),
                        thrust::raw_pointer_cast(_device_data->_d_force_ljcoul_z.data()));

	auto end_rbl_force = std::chrono::high_resolution_clock::now();
	std::chrono::duration<rbmd::Real> duration_rbl_force = end_rbl_force - start_rbl_force;
	TimingStatistics::Instance().record("Short-Range",duration_rbl_force.count());
  
    //energy
   if ("yes" == _energy_rbl_flag){
       ComputeLJCoulEnergy();
    }
}

void LJCutCoulKspace::ComputeLJVerlet()
{
  //neighbor_list_build
  auto start = std::chrono::high_resolution_clock::now();
  _list = _neighbor_list_builder->Build();

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration = end - start;
  TimingStatistics::Instance().record("Neighbor-List",duration.count());
  
  //
  auto start_verlet_force = std::chrono::high_resolution_clock::now();
  
  thrust::device_vector<rbmd::Real> d_total_evdwl(1, 0.0);
  thrust::device_vector<rbmd::Real> d_total_ecoul(1, 0.0);
  //
  auto num_atoms = *(_structure_info_data->_num_atoms);
  op::LJCutCoulForceOp<device::DEVICE_GPU>()(
                    *_box,_device_data->_d_erf_table, _cut_off, num_atoms,_alpha,_qqr2e,
                    thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
                    thrust::raw_pointer_cast(_device_data->_d_sigma.data()),
                    thrust::raw_pointer_cast(_device_data->_d_eps.data()),
                    thrust::raw_pointer_cast(_list->_start_idx.data()),
                    thrust::raw_pointer_cast(_list->_end_idx.data()),
                    thrust::raw_pointer_cast(_list->_d_neighbors.data()),
                    thrust::raw_pointer_cast(_device_data->_d_charge.data()),
                    thrust::raw_pointer_cast(_device_data->_d_px.data()),
                    thrust::raw_pointer_cast(_device_data->_d_py.data()),
                    thrust::raw_pointer_cast(_device_data->_d_pz.data()),
                    thrust::raw_pointer_cast(_device_data->_d_force_ljcoul_x.data()),
                    thrust::raw_pointer_cast(_device_data->_d_force_ljcoul_y.data()),
                    thrust::raw_pointer_cast(_device_data->_d_force_ljcoul_z.data()),
                    thrust::raw_pointer_cast(_device_data->_d_flat_virial_lj.data()),
                    thrust::raw_pointer_cast(d_total_evdwl.data()),
                      thrust::raw_pointer_cast(d_total_ecoul.data()));

  auto end_verlet_force = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration_verlet_force = end_verlet_force - start_verlet_force;
  TimingStatistics::Instance().record("Short-Range",duration_verlet_force.count());
  
  // 
  thrust::host_vector<rbmd::Real> h_total_evdwl(d_total_evdwl);
  thrust::host_vector<rbmd::Real> h_total_ecoul(d_total_ecoul);
  _e_vdwl = h_total_evdwl[0]/num_atoms;
  _e_coul = h_total_ecoul[0]/num_atoms;

  //sum virial_lj on host
  ReduceVirial(num_atoms,_device_data->_d_flat_virial_lj,
_device_data->_d_virial_lj);
}


void LJCutCoulKspace::SumForces()
{
  TransformForces(_device_data->_d_fx,_device_data->_d_force_ljcoul_x,
    _device_data->_d_force_kspace_x);

  TransformForces(_device_data->_d_fy,_device_data->_d_force_ljcoul_y,
    _device_data->_d_force_kspace_y);

  TransformForces(_device_data->_d_fz,_device_data->_d_force_ljcoul_z,
    _device_data->_d_force_kspace_z);
}


void LJCutCoulKspace::ComputeLJCoulEnergy()
{
  // energy
  //neighbor_list_build
  _list = _neighbor_list_builder->Build();

  //
  thrust::device_vector<rbmd::Real> _d_total_evdwl(1, 0.0);
  thrust::device_vector<rbmd::Real> _d_total_ecoul(1, 0.0);
  auto num_atoms = *(_structure_info_data->_num_atoms);
  op::LJCutCoulEnergyOp<device::DEVICE_GPU>()(
                *_box,_device_data->_d_erf_table,_cut_off,num_atoms,_alpha,_qqr2e,
                thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
                thrust::raw_pointer_cast(_device_data->_d_sigma.data()),
                thrust::raw_pointer_cast(_device_data->_d_eps.data()),
                thrust::raw_pointer_cast(_list->_start_idx.data()),
                thrust::raw_pointer_cast(_list->_end_idx.data()),
                thrust::raw_pointer_cast(_list->_d_neighbors.data()),
                thrust::raw_pointer_cast(_device_data->_d_charge.data()),
                thrust::raw_pointer_cast(_device_data->_d_px.data()),
                thrust::raw_pointer_cast(_device_data->_d_py.data()),
                thrust::raw_pointer_cast(_device_data->_d_pz.data()),
                thrust::raw_pointer_cast(_device_data->_d_flat_virial_lj.data()),
                thrust::raw_pointer_cast(_d_total_evdwl.data()),
                thrust::raw_pointer_cast(_d_total_ecoul.data()));

  // D2H
  thrust::host_vector<rbmd::Real> h_total_evdwl(_d_total_evdwl);
  thrust::host_vector<rbmd::Real> h_total_ecoul(_d_total_ecoul);
  _e_vdwl = h_total_evdwl[0]/num_atoms;
  _e_coul = h_total_ecoul[0]/num_atoms;

  //sum virial_lj on host
  ReduceVirial(num_atoms,_device_data->_d_flat_virial_lj,
_device_data->_d_virial_lj);
}

void LJCutCoulKspace::EvaluatePotentialenergy()
{
  _e_pe_rbl = _e_vdwl_rbl + _e_coul_rbl +_e_kspace;
  //test_ave_pe_rbl = _ave_pe_rbl;

  _e_pe = _e_vdwl+ _e_coul +_e_kspace;
  //test_ave_pe = _ave_pe;

  ThermoStats::Instance().AddThermoData("total-potential-energy",_e_pe);
  
  //out
  auto interval = DataManager::getInstance().getConfigData()->Get<rbmd::Id>(
"interval", "outputs", "thermo_out");

  std::ofstream outfile("thermo.txt", std::ios::app);
  if (outfile.tellp() == 0) {
    outfile << "step e_vdwl e_coul e_kspace e_pe" << std::endl;
  }
  if (test_current_step % interval == 0) {
    outfile << test_current_step << " " << _e_vdwl  << " "<< _e_coul <<" "
      << _e_kspace  << " " << _e_pe<< std::endl;
  }

  outfile.close();
}