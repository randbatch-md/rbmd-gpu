#include "cvff.h"

//#include <output/include/Logger.hpp>

#include "../../common/device_types.h"
#include "../../common/rbmd_define.h"
#include "../../common/types.h"
#include "../common/RBEPSample.h"
#include "../common/unit_factor.h"
#include "force_op/force_op.h"
#include "cvff_op/cvff_op.h"

#include "lj_op/lj_op.h"
#include "neighbor_list/include/linked_cell/linked_cell_locator.h"
#include "neighbor_list/include/neighbor_list_builder/full_neighbor_list_builder.h"
#include "neighbor_list/include/neighbor_list_builder/rbl_full_neighbor_list_builder.h"
#include "thrust/sort.h"
// #include <hipcub/hipcub.hpp>
// #include <hipcub/backend/rocprim/block/block_reduce.hpp>
#include "common/thermo_stats.hpp"
#include "common/timing_statistics.hpp"
extern int test_current_step;
extern std::map<std::string, UNIT> unit_factor_map;

CVFF::CVFF()
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

CVFF::~CVFF(){}

void CVFF::Init()
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

void CVFF::Execute() {
  const auto& config = DataManager::getInstance().getConfigData();

  ComputeLJCutCoulForce();
  if(config->PathExists({"hyper_parameters", "coulomb"})) {
    // 2. 委托 K-Space 组件计算长程力
    _kspace_calculator->Execute();
  }

  ComputeBondForce();
  ComputeAngleForce();
  if(*(_structure_info_data->_num_dihedrals)) {
    ComputeDihedralForce();
  }
  if(*(_structure_info_data->_num_impropers)) {
    ComputeImproperForce();
  }

  SumForces();

  EvaluatePotentialenergy();
}

void CVFF::ComputeLJCutCoulForce()
{
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

void CVFF::ComputeLJRBL()
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
    op::SpecialLJCutCoulRBLForceOp<device::DEVICE_GPU>()(
        *_box,_device_data->_d_erf_table,
      r_core, _cut_off, num_atoms,neighbor_sample_num,
      _rbl_list->_selection_frequency,_alpha,_qqr2e,
        thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
        thrust::raw_pointer_cast(_device_data->_d_atoms_id.data()),
        thrust::raw_pointer_cast(_device_data->_d_sigma.data()),
        thrust::raw_pointer_cast(_device_data->_d_eps.data()),
        thrust::raw_pointer_cast(_rbl_list->_start_idx.data()),
        thrust::raw_pointer_cast(_rbl_list->_end_idx.data()),
        thrust::raw_pointer_cast(_rbl_list->_d_neighbors.data()),
        thrust::raw_pointer_cast(_rbl_list->_d_random_neighbor.data()),
        thrust::raw_pointer_cast(_rbl_list->_d_random_neighbor_num.data()),
        thrust::raw_pointer_cast(_device_data->_d_special_ids.data()),
        thrust::raw_pointer_cast(_device_data->_d_special_weights.data()),
        thrust::raw_pointer_cast(_device_data->_d_special_offsets.data()),
        thrust::raw_pointer_cast(_device_data->_d_special_count.data()),
        thrust::raw_pointer_cast(_device_data->_d_charge.data()),
        thrust::raw_pointer_cast(_device_data->_d_px.data()),
        thrust::raw_pointer_cast(_device_data->_d_py.data()),
        thrust::raw_pointer_cast(_device_data->_d_pz.data()),
        thrust::raw_pointer_cast(_device_data->_d_force_ljcoul_x.data()),
        thrust::raw_pointer_cast(_device_data->_d_force_ljcoul_y.data()),
        thrust::raw_pointer_cast(_device_data->_d_force_ljcoul_z.data()));

    _corr_value_x =
        thrust::reduce(_device_data->_d_force_ljcoul_x.begin(), _device_data->_d_force_ljcoul_x.end(),
                       0.0f, thrust::plus<rbmd::Real>()) /num_atoms;
    _corr_value_y =
        thrust::reduce(_device_data->_d_force_ljcoul_y.begin(), _device_data->_d_force_ljcoul_y.end(),
                       0.0f, thrust::plus<rbmd::Real>()) /num_atoms;
    _corr_value_z =
        thrust::reduce(_device_data->_d_force_ljcoul_z.begin(), _device_data->_d_force_ljcoul_z.end(),
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
    if ("yes" == _energy_rbl_flag ) {
      ComputeLJCoulEnergy();
    }
}

void CVFF::ComputeLJVerlet()
{
  //neighbor_list_build
  auto start = std::chrono::high_resolution_clock::now();
  _list = _neighbor_list_builder->Build();

  auto end = std::chrono::high_resolution_clock::now();

  std::chrono::duration<rbmd::Real> duration = end - start;

  TimingStatistics::Instance().record("Neighbor-List",duration.count());

  //
  auto start_verlet_force = std::chrono::high_resolution_clock::now();
  //
  thrust::device_vector<rbmd::Real> _d_total_evdwl(1, 0.0);
  thrust::device_vector<rbmd::Real> _d_total_ecoul(1, 0.0);
  //
  auto num_atoms = *(_structure_info_data->_num_atoms);
  op::SpecialLJCutCoulForceOp<device::DEVICE_GPU>()(
                  *_box,_device_data->_d_erf_table, _cut_off, num_atoms,_alpha,_qqr2e,
                  thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
                  thrust::raw_pointer_cast(_device_data->_d_atoms_id.data()),
                  thrust::raw_pointer_cast(_device_data->_d_sigma.data()),
                  thrust::raw_pointer_cast(_device_data->_d_eps.data()),
                  thrust::raw_pointer_cast(_list->_start_idx.data()),
                  thrust::raw_pointer_cast(_list->_end_idx.data()),
                  thrust::raw_pointer_cast(_list->_d_neighbors.data()),
                  thrust::raw_pointer_cast(_device_data->_d_special_ids.data()),
                  thrust::raw_pointer_cast(_device_data->_d_special_weights.data()),
                  thrust::raw_pointer_cast(_device_data->_d_special_offsets.data()),
                  thrust::raw_pointer_cast(_device_data->_d_special_count.data()),
                  thrust::raw_pointer_cast(_device_data->_d_charge.data()),
                  thrust::raw_pointer_cast(_device_data->_d_px.data()),
                  thrust::raw_pointer_cast(_device_data->_d_py.data()),
                  thrust::raw_pointer_cast(_device_data->_d_pz.data()),
                  thrust::raw_pointer_cast(_device_data->_d_force_ljcoul_x.data()),
                  thrust::raw_pointer_cast(_device_data->_d_force_ljcoul_y.data()),
                  thrust::raw_pointer_cast(_device_data->_d_force_ljcoul_z.data()),
                  thrust::raw_pointer_cast(_device_data->_d_flat_virial_lj.data()),
                  thrust::raw_pointer_cast(_d_total_evdwl.data()),
                  thrust::raw_pointer_cast(_d_total_ecoul.data()));

  auto end_verlet_force = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration_verlet_force = end_verlet_force - start_verlet_force;
  TimingStatistics::Instance().record("Short-Range",duration_verlet_force.count());
  // D2H
  thrust::host_vector<rbmd::Real> h_total_evdwl(_d_total_evdwl);
  thrust::host_vector<rbmd::Real> h_total_ecoul(_d_total_ecoul);
  _e_vdwl = h_total_evdwl[0]/num_atoms;
  _e_coul = h_total_ecoul[0]/num_atoms;

//sum virial_special_lj on host
  ReduceVirial(num_atoms,_device_data->_d_flat_virial_lj,
_device_data->_d_virial_lj);
}

void CVFF::ComputeLJCoulEnergy()
{
  // energy
  //neighbor_list_build
  //auto start = std::chrono::high_resolution_clock::now();
  _list = _neighbor_list_builder->Build();

  //auto end = std::chrono::high_resolution_clock::now();

  //std::chrono::duration<rbmd::Real> duration = end - start;
  //TimingStatistics::Instance().record("Neighbor-List",duration.count());


  thrust::device_vector<rbmd::Real> _d_total_evdwl(1, 0.0);
  thrust::device_vector<rbmd::Real> _d_total_ecoul(1, 0.0);

  auto num_atoms = *(_structure_info_data->_num_atoms);
  op::SpeciaLJCutCoulEnergyOp<device::DEVICE_GPU>()(
                *_box,_device_data->_d_erf_table,_cut_off, num_atoms,_alpha,_qqr2e,
                thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
                thrust::raw_pointer_cast(_device_data->_d_atoms_id.data()),
                thrust::raw_pointer_cast(_device_data->_d_sigma.data()),
                thrust::raw_pointer_cast(_device_data->_d_eps.data()),
                thrust::raw_pointer_cast(_list->_start_idx.data()),
                thrust::raw_pointer_cast(_list->_end_idx.data()),
                thrust::raw_pointer_cast(_list->_d_neighbors.data()),
                thrust::raw_pointer_cast(_device_data->_d_special_ids.data()),
                thrust::raw_pointer_cast(_device_data->_d_special_weights.data()),
                thrust::raw_pointer_cast(_device_data->_d_special_offsets.data()),
                thrust::raw_pointer_cast(_device_data->_d_special_count.data()),
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

  //sum virial on host
  ReduceVirial(num_atoms,_device_data->_d_flat_virial_lj,
_device_data->_d_virial_lj);
}

void CVFF::SumForces()
{
  TransformForces(_device_data->_d_fx,_device_data->_d_force_ljcoul_x,
    _device_data->_d_force_kspace_x,_device_data->_d_force_bond_x,
    _device_data->_d_force_angle_x,_device_data->_d_force_dihedral_x,
    _device_data->_d_force_improper_x);

  TransformForces(_device_data->_d_fy,_device_data->_d_force_ljcoul_y,
    _device_data->_d_force_kspace_y,_device_data->_d_force_bond_y,
    _device_data->_d_force_angle_y,_device_data->_d_force_dihedral_y,
    _device_data->_d_force_improper_y);

  TransformForces(_device_data->_d_fz,_device_data->_d_force_ljcoul_z,
    _device_data->_d_force_kspace_z,_device_data->_d_force_bond_z,
    _device_data->_d_force_angle_z,_device_data->_d_force_dihedral_z,
    _device_data->_d_force_improper_z);
}

void CVFF::ComputeBondForce()
{
  auto start = std::chrono::high_resolution_clock::now();

  auto _atom_id_to_idx =
    LinkedCellLocator::GetInstance().GetLinkedCell()->_atom_id_to_idx;

  thrust::fill(_device_data->_d_force_bond_x.begin(),
    _device_data->_d_force_bond_x.end(), 0.0f);
  thrust::fill(_device_data->_d_force_bond_y.begin(),
    _device_data->_d_force_bond_y.end(), 0.0f);
  thrust::fill(_device_data->_d_force_bond_z.begin(),
    _device_data->_d_force_bond_z.end(), 0.0f);

  thrust::fill(_device_data->_d_flat_virial_bond_atom.begin(),
    _device_data->_d_flat_virial_bond_atom.end(), 0.0f);

  thrust::device_vector<rbmd::Real> d_total_ebond(1, 0.0);

  auto num_bonds = *(_structure_info_data->_num_bonds);
  auto num_atoms = *(_structure_info_data->_num_atoms);
  op::ComputeBondForceOp<device::DEVICE_GPU>()(
    *_box,num_atoms,num_bonds,thrust::raw_pointer_cast(_atom_id_to_idx.data()),
    thrust::raw_pointer_cast(_device_data->_d_bond_coeffs_k.data()),
    thrust::raw_pointer_cast(_device_data->_d_bond_coeffs_equilibrium.data()),
    thrust::raw_pointer_cast(_device_data->_d_bond_type.data()),
    thrust::raw_pointer_cast(_device_data->_d_bond_id0.data()),
    thrust::raw_pointer_cast(_device_data->_d_bond_id1.data()),
    thrust::raw_pointer_cast(_device_data->_d_px.data()),
    thrust::raw_pointer_cast(_device_data->_d_py.data()),
    thrust::raw_pointer_cast(_device_data->_d_pz.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_bond_x.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_bond_y.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_bond_z.data()),
    thrust::raw_pointer_cast(_device_data->_d_flat_virial_bond_atom.data()),
    thrust::raw_pointer_cast(_device_data->_d_flat_virial_bond_list.data()),
    thrust::raw_pointer_cast(d_total_ebond.data()));

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration = end - start;
  TimingStatistics::Instance().record("Bond",duration.count());

  // D2H
  thrust::host_vector<rbmd::Real> h_total_ebond(d_total_ebond);
  _e_bond = h_total_ebond[0]/num_bonds;

  ThermoStats::Instance().AddThermoData("bond",_e_bond);

  // //sum virial_bond  on host
  ReduceVirial(num_atoms,_device_data->_d_flat_virial_bond_atom,
_device_data->_d_virial_bond);
}

void CVFF::ComputeAngleForce()
{
  auto start = std::chrono::high_resolution_clock::now();
  auto atom_id_to_idx =
    LinkedCellLocator::GetInstance().GetLinkedCell()->_atom_id_to_idx;

  thrust::fill(_device_data->_d_force_angle_x.begin(),
  _device_data->_d_force_angle_x.end(), 0.0f);
  thrust::fill(_device_data->_d_force_angle_y.begin(),
    _device_data->_d_force_angle_y.end(), 0.0f);
  thrust::fill(_device_data->_d_force_angle_z.begin(),
    _device_data->_d_force_angle_z.end(), 0.0f);

  thrust::fill(_device_data->_d_flat_virial_angle_atom.begin(),
  _device_data->_d_flat_virial_angle_atom.end(), 0.0f);

  thrust::device_vector<rbmd::Real> d_total_eangle(1, 0.0);

  auto num_angles = *(_structure_info_data->_num_angles);
  auto num_atoms = *(_structure_info_data->_num_atoms);
  op::ComputeAngleForceOp<device::DEVICE_GPU>()(
    *_box,num_atoms,num_angles,
    thrust::raw_pointer_cast(atom_id_to_idx.data()),
    thrust::raw_pointer_cast(_device_data->_d_angle_coeffs_k.data()),
    thrust::raw_pointer_cast(_device_data->_d_angle_coeffs_equilibrium.data()),
    thrust::raw_pointer_cast(_device_data->_d_angle_type.data()),
    thrust::raw_pointer_cast(_device_data->_d_angle_id0.data()),
    thrust::raw_pointer_cast(_device_data->_d_angle_id1.data()),
    thrust::raw_pointer_cast(_device_data->_d_angle_id2.data()),
    thrust::raw_pointer_cast(_device_data->_d_px.data()),
    thrust::raw_pointer_cast(_device_data->_d_py.data()),
    thrust::raw_pointer_cast(_device_data->_d_pz.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_angle_x.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_angle_y.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_angle_z.data()),
    thrust::raw_pointer_cast(_device_data->_d_flat_virial_angle_atom.data()),
    thrust::raw_pointer_cast(_device_data->_d_flat_virial_angle_list.data()),
    thrust::raw_pointer_cast(d_total_eangle.data()));

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration = end - start;
  TimingStatistics::Instance().record("Angle",duration.count());
  // D2H
  thrust::host_vector<rbmd::Real> h_total_eangle(d_total_eangle);
  _e_angle = h_total_eangle[0]/num_angles;

  ThermoStats::Instance().AddThermoData("angle",_e_angle);

  //sum virial_angle on host
  ReduceVirial(num_atoms,_device_data->_d_flat_virial_angle_atom,
_device_data->_d_virial_angle);
}

void CVFF::ComputeDihedralForce()
{
  auto dihedral_type = DataManager::getInstance().getConfigData()->Get
  <std::string>("dihedral_type", "hyper_parameters", "force_field");

  if (dihedral_type == "harmonic") {
    DihedralHarmonic();
  }
  else if (dihedral_type == "opls") {
    DihedralOPLS();
  }

  //add thermo
  ThermoStats::Instance().AddThermoData("dihedral",_e_dihedral);
}

void CVFF::DihedralHarmonic() {
   auto start = std::chrono::high_resolution_clock::now();
  thrust::fill(_device_data->_d_force_dihedral_x.begin(),
    _device_data->_d_force_dihedral_x.end(), 0.0f);
  thrust::fill(_device_data->_d_force_dihedral_y.begin(),
    _device_data->_d_force_dihedral_y.end(), 0.0f);
  thrust::fill(_device_data->_d_force_dihedral_z.begin(),
    _device_data->_d_force_dihedral_z.end(), 0.0f);

  thrust::fill(_device_data->_d_flat_virial_dihedral_atom.begin(),
  _device_data->_d_flat_virial_dihedral_atom.end(), 0.0f);

  //thrust::device_vector<int4> dihedral_list;
  auto atom_id_to_idx =
    LinkedCellLocator::GetInstance().GetLinkedCell()->_atom_id_to_idx;

  thrust::device_vector<rbmd::Real> d_total_edihedral(1, 0.0);

  auto num_atoms = *(_structure_info_data->_num_atoms);
  auto num_dihedrals = *(_structure_info_data->_num_dihedrals);
  op::ComputeDihedralForceOp<device::DEVICE_GPU>()(
    *_box,num_atoms,num_dihedrals,
    thrust::raw_pointer_cast(atom_id_to_idx.data()),
    thrust::raw_pointer_cast(_device_data->_d_dihedral_coeffs_k.data()),
    thrust::raw_pointer_cast(_device_data->_d_dihedral_coeffs_sign.data()),
    thrust::raw_pointer_cast(_device_data->_d_dihedral_coeffs_multiplicity.data()),
    thrust::raw_pointer_cast(_device_data->_d_dihedral_type.data()),
    thrust::raw_pointer_cast(_device_data->_d_dihedral_id0.data()),
    thrust::raw_pointer_cast(_device_data->_d_dihedral_id1.data()),
    thrust::raw_pointer_cast(_device_data->_d_dihedral_id2.data()),
    thrust::raw_pointer_cast(_device_data->_d_dihedral_id3.data()),
    thrust::raw_pointer_cast(_device_data->_d_px.data()),
    thrust::raw_pointer_cast(_device_data->_d_py.data()),
    thrust::raw_pointer_cast(_device_data->_d_pz.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_dihedral_x.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_dihedral_y.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_dihedral_z.data()),
    thrust::raw_pointer_cast(_device_data->_d_flat_virial_dihedral_atom.data()),
    thrust::raw_pointer_cast(_device_data->_d_flat_virial_dihedral_list.data()),
    thrust::raw_pointer_cast(d_total_edihedral.data()));

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration = end - start;
  TimingStatistics::Instance().record("Dihedral",duration.count());

  // D2H
  thrust::host_vector<rbmd::Real> h_total_edihedral(d_total_edihedral);
  _e_dihedral = h_total_edihedral[0]/num_dihedrals;

  //sum virial_dihedral on host
  ReduceVirial(num_atoms,_device_data->_d_flat_virial_dihedral_atom,
  _device_data->_d_virial_dihedral);
}

void CVFF::DihedralOPLS() {
   auto start = std::chrono::high_resolution_clock::now();
  thrust::fill(_device_data->_d_force_dihedral_x.begin(),
    _device_data->_d_force_dihedral_x.end(), 0.0f);
  thrust::fill(_device_data->_d_force_dihedral_y.begin(),
    _device_data->_d_force_dihedral_y.end(), 0.0f);
  thrust::fill(_device_data->_d_force_dihedral_z.begin(),
    _device_data->_d_force_dihedral_z.end(), 0.0f);

  thrust::fill(_device_data->_d_flat_virial_dihedral_atom.begin(),
  _device_data->_d_flat_virial_dihedral_atom.end(), 0.0f);

  //thrust::device_vector<int4> dihedral_list;
  auto atom_id_to_idx =
    LinkedCellLocator::GetInstance().GetLinkedCell()->_atom_id_to_idx;

  thrust::device_vector<rbmd::Real> d_total_edihedral(1, 0.0);

  auto num_atoms = *(_structure_info_data->_num_atoms);
  auto num_dihedrals = *(_structure_info_data->_num_dihedrals);
  op::ComputeDihedralOPLSForceOp<device::DEVICE_GPU>()(
    *_box,num_atoms,num_dihedrals,
    thrust::raw_pointer_cast(atom_id_to_idx.data()),
    thrust::raw_pointer_cast(_device_data->_d_dihedral_coeffs_k1.data()),
    thrust::raw_pointer_cast(_device_data->_d_dihedral_coeffs_k2.data()),
    thrust::raw_pointer_cast(_device_data->_d_dihedral_coeffs_k3.data()),
    thrust::raw_pointer_cast(_device_data->_d_dihedral_coeffs_k4.data()),
    thrust::raw_pointer_cast(_device_data->_d_dihedral_type.data()),
    thrust::raw_pointer_cast(_device_data->_d_dihedral_id0.data()),
    thrust::raw_pointer_cast(_device_data->_d_dihedral_id1.data()),
    thrust::raw_pointer_cast(_device_data->_d_dihedral_id2.data()),
    thrust::raw_pointer_cast(_device_data->_d_dihedral_id3.data()),
    thrust::raw_pointer_cast(_device_data->_d_px.data()),
    thrust::raw_pointer_cast(_device_data->_d_py.data()),
    thrust::raw_pointer_cast(_device_data->_d_pz.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_dihedral_x.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_dihedral_y.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_dihedral_z.data()),
    thrust::raw_pointer_cast(_device_data->_d_flat_virial_dihedral_atom.data()),
    thrust::raw_pointer_cast(_device_data->_d_flat_virial_dihedral_list.data()),
    thrust::raw_pointer_cast(d_total_edihedral.data()));

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration = end - start;
  TimingStatistics::Instance().record("Dihedral",duration.count());
  // D2H
  thrust::host_vector<rbmd::Real> h_total_edihedral(d_total_edihedral);
  _e_dihedral = h_total_edihedral[0]/num_dihedrals;

  //sum virial_dihedral on host
  ReduceVirial(num_atoms,_device_data->_d_flat_virial_dihedral_atom,
    _device_data->_d_virial_dihedral);
}

void CVFF::ComputeImproperForce()
{
  auto improper_type = DataManager::getInstance().getConfigData()->Get
<std::string>("improper_type", "hyper_parameters", "force_field");

  if (improper_type == "harmonic") {
    ImproperHarmonic();
  }
  else if (improper_type == "cvff") {
    ImproperCVFF();
  }

  //add thermo
  ThermoStats::Instance().AddThermoData("improper",_e_improper);
}

void CVFF::ImproperHarmonic() {
  auto start = std::chrono::high_resolution_clock::now();
  thrust::fill(_device_data->_d_force_improper_x.begin(),
    _device_data->_d_force_improper_x.end(), 0.0f);
  thrust::fill(_device_data->_d_force_improper_y.begin(),
    _device_data->_d_force_improper_y.end(), 0.0f);
  thrust::fill(_device_data->_d_force_improper_z.begin(),
    _device_data->_d_force_improper_z.end(), 0.0f);

  auto atom_id_to_idx =
    LinkedCellLocator::GetInstance().GetLinkedCell()->_atom_id_to_idx;

  thrust::device_vector<rbmd::Real> d_total_eimproper(1, 0.0);

  auto num_atoms = *(_structure_info_data->_num_atoms);
  auto num_impropers = *(_structure_info_data->_num_impropers);
  op::ComputeImproperHarmonicForceOp<device::DEVICE_GPU>()(
    *_box,num_atoms,num_impropers,
    thrust::raw_pointer_cast(atom_id_to_idx.data()),
    thrust::raw_pointer_cast(_device_data->_d_improper_coeffs_k.data()),
    thrust::raw_pointer_cast(_device_data->_d_improper_coeffs_chi.data()),
    thrust::raw_pointer_cast(_device_data->_d_improper_type.data()),
    thrust::raw_pointer_cast(_device_data->_d_improper_id0.data()),
    thrust::raw_pointer_cast(_device_data->_d_improper_id1.data()),
    thrust::raw_pointer_cast(_device_data->_d_improper_id2.data()),
    thrust::raw_pointer_cast(_device_data->_d_improper_id3.data()),
    thrust::raw_pointer_cast(_device_data->_d_px.data()),
    thrust::raw_pointer_cast(_device_data->_d_py.data()),
    thrust::raw_pointer_cast(_device_data->_d_pz.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_improper_x.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_improper_y.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_improper_z.data()),
    thrust::raw_pointer_cast(_device_data->_d_flat_virial_improper_atom.data()),
    thrust::raw_pointer_cast(d_total_eimproper.data()));

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration = end - start;
  TimingStatistics::Instance().record("Improper",duration.count());

  // D2H
  thrust::host_vector<rbmd::Real> h_total_eimproper(d_total_eimproper);
  _e_improper = h_total_eimproper[0]/num_impropers;

  ReduceVirial(num_atoms,_device_data->_d_flat_virial_improper_atom,
  _device_data->_d_virial_improper);
}

void CVFF::ImproperCVFF()
{
   auto start = std::chrono::high_resolution_clock::now();
  thrust::fill(_device_data->_d_force_improper_x.begin(),
    _device_data->_d_force_improper_x.end(), 0.0f);
  thrust::fill(_device_data->_d_force_improper_y.begin(),
    _device_data->_d_force_improper_y.end(), 0.0f);
  thrust::fill(_device_data->_d_force_improper_z.begin(),
    _device_data->_d_force_improper_z.end(), 0.0f);

  auto atom_id_to_idx =
    LinkedCellLocator::GetInstance().GetLinkedCell()->_atom_id_to_idx;

  thrust::device_vector<rbmd::Real> d_total_eimproper(1, 0.0);

  auto num_atoms = *(_structure_info_data->_num_atoms);
  auto num_impropers = *(_structure_info_data->_num_impropers);
  op::ComputeImproperCVFFForceOp<device::DEVICE_GPU>()(
    *_box,num_atoms,num_impropers,
    thrust::raw_pointer_cast(atom_id_to_idx.data()),
    thrust::raw_pointer_cast(_device_data->_d_improper_coeffs_k.data()),
    thrust::raw_pointer_cast(_device_data->_d_improper_coeffs_d.data()),
    thrust::raw_pointer_cast(_device_data->_d_improper_coeffs_n.data()),
    thrust::raw_pointer_cast(_device_data->_d_improper_type.data()),
    thrust::raw_pointer_cast(_device_data->_d_improper_id0.data()),
    thrust::raw_pointer_cast(_device_data->_d_improper_id1.data()),
    thrust::raw_pointer_cast(_device_data->_d_improper_id2.data()),
    thrust::raw_pointer_cast(_device_data->_d_improper_id3.data()),
    thrust::raw_pointer_cast(_device_data->_d_px.data()),
    thrust::raw_pointer_cast(_device_data->_d_py.data()),
    thrust::raw_pointer_cast(_device_data->_d_pz.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_improper_x.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_improper_y.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_improper_z.data()),
    thrust::raw_pointer_cast(_device_data->_d_flat_virial_improper_atom.data()),
    thrust::raw_pointer_cast(d_total_eimproper.data()));

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration = end - start;
  TimingStatistics::Instance().record("Improper",duration.count());

  // D2H
  thrust::host_vector<rbmd::Real> h_total_eimproper(d_total_eimproper);
  _e_improper = h_total_eimproper[0]/num_impropers;

  ReduceVirial(num_atoms,_device_data->_d_flat_virial_improper_atom,
_device_data->_d_virial_improper);
}

void CVFF::EvaluatePotentialenergy()
{
  _e_pe_rbl = _e_vdwl_rbl + _e_coul_rbl+_e_kspace+
                  _e_bond + _e_angle+_e_dihedral+_e_improper;

  _e_pe = _e_vdwl+ _e_coul +_e_kspace+
              _e_bond +_e_angle +_e_dihedral+_e_improper;

  ThermoStats::Instance().AddThermoData("total-potential-energy",_e_pe);

  //out
  auto interval = DataManager::getInstance().getConfigData()->Get<rbmd::Id>(
"interval", "outputs", "thermo_out");

  std::ofstream outfile("thermo.txt", std::ios::app);
  if (outfile.tellp() == 0) {
    outfile << "step  e_vdwl  e_coul  e_kspace  e_bond "
            << "e_angle  e_dihedral e_improper e_pe" << std::endl;
  }
  if (test_current_step % interval == 0) {
    outfile << test_current_step << " " << _e_vdwl << " " << _e_coul << " "
        << _e_kspace << " " << _e_bond << " " << _e_angle << " "
        << _e_dihedral << " " << _e_improper << " " <<  _e_pe << std::endl;
  }
  outfile.close();
}



