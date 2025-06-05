#include "lj_cut_coul_kspace.h"

#include "../../common/device_types.h"
#include "../../common/rbmd_define.h"
#include "../../common/types.h"
#include "../../common/unit_factor.h"
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

  auto unit = DataManager::getInstance().getConfigData()->Get
<std::string>("unit", "init_configuration", "read_data");
  UNIT unit_factor = unit_factor_map[unit];
  switch (unit_factor) {
    case UNIT::LJ:
      _qqr2e = UnitFactor<UNIT::LJ>::_qqr2e;
    break;

    case UNIT::REAL:
      _qqr2e = UnitFactor<UNIT::REAL>::_qqr2e;
    break;

    default:
      break;
  }

//   _accuracy = DataManager::getInstance().getConfigData()->Get<rbmd::Real>(
// "accuracy", "hyper_parameters", "coulomb");
//   auto box =  DataManager::getInstance().getMDData()->_box;
//   auto volue = CalculateVolume(*box);
//   auto num_atoms = *(_structure_info_data->_num_atoms);
//   //  sum q_sq
//   ComputeQsqSum(); //q2
//
//   //compute g_ewald
//   _g_ewald = _accuracy*SQRT(num_atoms*_cut_off*volue) / (2.0*_sum_sq_charge);
//   if (_g_ewald >= 1.0) _g_ewald = (1.35 - 0.15*LOG(_accuracy))/_cut_off;
//   else _g_ewald = SQRT(-LOG(_g_ewald)) / _cut_off;
//   _alpha = _g_ewald*_g_ewald;
//
//   //automatically compute kmax
//   SetKspacePara(); //kmax
//   //ComputeEwlad_fix();
//
//   std::cout << "g_ewald: " << _g_ewald  <<", alpha: "<<
//     _alpha  << ", num_k: " << _num_k << std::endl;
  std::remove("thermo.txt");
}

LJCutCoulKspace::~LJCutCoulKspace(){}

void LJCutCoulKspace::Init()
{
  _cut_off = DataManager::getInstance().getConfigData()->Get
      <rbmd::Real>("cut_off", "hyper_parameters", "neighbor");

   _neighbor_type = DataManager::getInstance().getConfigData()->Get
      <std::string>("type", "hyper_parameters", "neighbor");

   _coulomb_type =DataManager::getInstance().getConfigData()->Get<std::string>(
        "type", "hyper_parameters", "coulomb");

  _alpha = DataManager::getInstance().getConfigData()->Get<rbmd::Real>(
"alpha", "hyper_parameters", "coulomb");
  auto Kmax =DataManager::getInstance().getConfigData()->
    GetArray<rbmd::Id>("kmax", "hyper_parameters", "coulomb");
  _kmax_array.x = Kmax[0];
  _kmax_array.y = Kmax[1];
  _kmax_array.z = Kmax[2];
  _num_k =  (2*_kmax_array.x +1)  * (2*_kmax_array.y +1)
            * (2*_kmax_array.z +1) - 1;

  if("RBE" == _coulomb_type) {
    _RBE_P = DataManager::getInstance().getConfigData()->Get<rbmd::Id>(
  "coulomb_sample_num", "hyper_parameters", "coulomb");
    GetPsampleKey();

    _energy_rbe_flag = DataManager::getInstance().getConfigData()->Get<std::string>("energy_rbe_flag", "hyper_parameters", "coulomb");
  }
}

void LJCutCoulKspace::Execute()
{
  ComputeLJCutCoulForce();
  ComputeKspaceForce();
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
   _energy_rbl_flag = DataManager::getInstance().getConfigData()->Get<std::string>
      ("energy_rbl_flag", "hyper_parameters", "neighbor");
   if ("yes" == _energy_rbl_flag ){
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

void LJCutCoulKspace::ComputeKspaceForce()
{
  if("RBE" == _coulomb_type)
  {
      ComputeRBE();
  }
  else
  {
     ComputeEwlad();
  }
  //add thermo
  ThermoStats::Instance().AddThermoData("kspace",_e_kspace);
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

void LJCutCoulKspace::ComputeChargeStructureFactorEwald(
    Box box,
    rbmd::Id num_atoms,
    Int3 kmax_array,
    rbmd::Real alpha,
    rbmd::Real qqr2e,
    thrust::host_vector<rbmd::Real> value_Re_array,
    thrust::host_vector<rbmd::Real> value_Im_array)
{
    //thrust::fill(density_real.begin(), density_real.end(), 0.0f);
    //thrust::fill(density_imag.begin(), density_imag.end(), 0.0f);
    auto start = std::chrono::high_resolution_clock::now();
    thrust::device_vector<rbmd::Real> density_real_atom;
    thrust::device_vector<rbmd::Real> density_imag_atom;
    density_real_atom.resize(num_atoms);
    density_imag_atom.resize(num_atoms);

    rbmd::Real total_energy_kspace= 0;
    rbmd::Id index = 0;
    for (rbmd::Id i = -INT_DATA(kmax_array)[0]; i <= INT_DATA(kmax_array)[0]; i++)
    {
        for (rbmd::Id j = -INT_DATA(kmax_array)[1]; j <= INT_DATA(kmax_array)[1]; j++)
        {
            for (rbmd::Id k = -INT_DATA(kmax_array)[2]; k <= INT_DATA(kmax_array)[2]; k++)
            {
                if (!(i == 0 && j == 0 && k == 0))
                {
                    Real3 K = make_Real3(rbmd::Real(2 * M_PI * i / box._length[0]),
                                           rbmd::Real(2 * M_PI * j / box._length[1]),
                                           rbmd::Real(2 * M_PI * k / box._length[2]));
                    rbmd::Real Range_K = SQRT(K.x * K.x + K.y * K.y + K.z * K.z);
                    rbmd::Real Range_K2 = Range_K*Range_K;
                    rbmd::Real alpha_inv =  1 / alpha;

                    op::ComputeChargeStructureFactorOp<device::DEVICE_GPU>()(
                        num_atoms, K,
                        thrust::raw_pointer_cast(_device_data->_d_charge.data()),
                        thrust::raw_pointer_cast(_device_data->_d_px.data()),
                        thrust::raw_pointer_cast(_device_data->_d_py.data()),
                        thrust::raw_pointer_cast(_device_data->_d_pz.data()),
                        thrust::raw_pointer_cast(density_real_atom.data()),
                        thrust::raw_pointer_cast(density_imag_atom.data()));

                    rbmd::Real value_Re = thrust::reduce(density_real_atom.begin(),
                      density_real_atom.end(), 0.0f, thrust::plus<rbmd::Real>());
                    rbmd::Real value_Im = thrust::reduce(density_imag_atom.begin(),
                      density_imag_atom.end(), 0.0f, thrust::plus<rbmd::Real>());
                    rbmd::Real Range_density2 = POW(value_Re, 2.0) + POW(value_Im, 2.0);

                    total_energy_kspace +=
                      EXP(-0.25 * Range_K2 * alpha_inv) * Range_density2 / Range_K2;

                    value_Re_array[index] = value_Re;
                    value_Im_array[index] = value_Im;
                    index++;
                }
            }
        }
    }

  //energy

  //charge self energy//
  ComputeSelfEnergy(alpha,qqr2e,_e_self_energy);

  //compute Kspace energy//
  rbmd::Real volume = box._length[0] * box._length[1]*box._length[2];
  total_energy_kspace = qqr2e * (2 * M_PI / volume) * total_energy_kspace;
  _e_kspace = total_energy_kspace / num_atoms;

  _e_kspace = _e_kspace + _e_self_energy;
}

void LJCutCoulKspace::ComputeEwlad()
{
  auto start = std::chrono::high_resolution_clock::now();
  //
  auto num_atoms = *(_structure_info_data->_num_atoms);

  //compute charge structure factor
  thrust::host_vector<rbmd::Real> h_Re_array(_num_k);
  thrust::host_vector<rbmd::Real> h_Im_array(_num_k);
  ComputeChargeStructureFactorEwald(*_box, num_atoms, _kmax_array,
    _alpha,_qqr2e, h_Re_array,h_Im_array);

  thrust::device_vector<rbmd::Real> d_real_array(_num_k);
  thrust::device_vector<rbmd::Real> d_imag_array(_num_k);
  d_real_array = h_Re_array;
  d_imag_array = h_Im_array;

  //EwaldForce//
  op::ComputeEwaldForceOp<device::DEVICE_GPU>()(
        *_box,num_atoms, _kmax_array, _alpha,_qqr2e,
        thrust::raw_pointer_cast(d_real_array.data()),
          thrust::raw_pointer_cast(d_imag_array.data()),
        thrust::raw_pointer_cast(_device_data->_d_charge.data()),
        thrust::raw_pointer_cast(_device_data->_d_px.data()),
        thrust::raw_pointer_cast(_device_data->_d_py.data()),
        thrust::raw_pointer_cast(_device_data->_d_pz.data()),
        thrust::raw_pointer_cast(_device_data->_d_force_kspace_x.data()),
        thrust::raw_pointer_cast(_device_data->_d_force_kspace_y.data()),
        thrust::raw_pointer_cast(_device_data->_d_force_kspace_z.data()),
        thrust::raw_pointer_cast(_device_data->_d_flat_virial_kspace.data()));
		
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration = end - start;
  TimingStatistics::Instance().record("Long-Range",duration.count());
  
  //sum virial_kspace on host
  ReduceVirial(num_atoms,_device_data->_d_flat_virial_kspace,
_device_data->_d_virial_kspace);
}


void LJCutCoulKspace::RBEInit(Box box,rbmd::Real alpha,rbmd::Id RBE_P)
{
  auto num_atoms = *(_structure_info_data->_num_atoms);
  Real3 sigma = { rbmd::Real((SQRT(alpha / 2.0) * box._length[0]/M_PI)),
                  rbmd::Real((SQRT(alpha / 2.0) * box._length[1]/M_PI)),
                  rbmd::Real((SQRT(alpha / 2.0) * box._length[2]/M_PI))};
  auto random = true;
  RBEPSAMPLE rbe_presolve_psample = { alpha, box, RBE_P, random};

  thrust::host_vector<rbmd::Real> h_P_Sample_x(RBE_P);
  thrust::host_vector<rbmd::Real> h_P_Sample_y(RBE_P);
  thrust::host_vector<rbmd::Real> h_P_Sample_z(RBE_P);

  rbe_presolve_psample.Fetch_P_Sample(0.0, sigma,
    thrust::raw_pointer_cast(h_P_Sample_x.data()),
    thrust::raw_pointer_cast(h_P_Sample_y.data()),
    thrust::raw_pointer_cast(h_P_Sample_z.data()));

  _P_Sample_x = h_P_Sample_x;
  _P_Sample_y = h_P_Sample_y;
  _P_Sample_z = h_P_Sample_z;

  //index key
  _psample_key.resize(num_atoms * RBE_P);
  op::GenerateIndexArrayOp<device::DEVICE_GPU>()(
    num_atoms,RBE_P,
    thrust::raw_pointer_cast(_psample_key.data()));
}

void LJCutCoulKspace::GetPsampleKey()
{
  //index key
  auto num_atoms = *(_structure_info_data->_num_atoms);
  _psample_key.resize(num_atoms * _RBE_P);
  op::GenerateIndexArrayOp<device::DEVICE_GPU>()(num_atoms,_RBE_P,
    thrust::raw_pointer_cast(_psample_key.data()));
}

void LJCutCoulKspace::ComputeChargeStructureFactorRBE(
   Box box,
   rbmd::Id num_atoms,
   Int3 kmax_array,
   rbmd::Real alpha,
   rbmd::Id RBE_P,
   rbmd::Real qqr2e,
   thrust::device_vector<rbmd::Real> rhok_real_redue,
   thrust::device_vector<rbmd::Real> rhok_image_redue)
{
  //get P_Sample
  RBEInit(*_box,_alpha,_RBE_P);

  thrust::device_vector<rbmd::Real>  rhok_real_atom;
  thrust::device_vector<rbmd::Real>  rhok_image_atom;
  rhok_real_atom.resize(num_atoms* RBE_P);
  rhok_image_atom.resize(num_atoms* RBE_P);
  auto p_number= RBE_P;

  //Charge Structure Factor
  op::ComputePnumberChargeStructureFactorOp<device::DEVICE_GPU>()(
      box, num_atoms, p_number,
      thrust::raw_pointer_cast(_device_data->_d_charge.data()),
      thrust::raw_pointer_cast(_P_Sample_x.data()),
      raw_pointer_cast(_P_Sample_y.data()),
      raw_pointer_cast(_P_Sample_z.data()),
      thrust::raw_pointer_cast(_device_data->_d_px.data()),
      thrust::raw_pointer_cast(_device_data->_d_py.data()),
      thrust::raw_pointer_cast(_device_data->_d_pz.data()),
      thrust::raw_pointer_cast(rhok_real_atom.data()),
      thrust::raw_pointer_cast(rhok_image_atom.data()));

  // reduce_by_key for rhok
  thrust::device_vector<rbmd::Id>  psamplekey_out;
  psamplekey_out.resize(RBE_P);

   reduce_by_key(_psample_key.begin(), _psample_key.end(),
    rhok_real_atom.begin(),psamplekey_out.begin(), rhok_real_redue.begin(),
    thrust::equal_to<rbmd::Id>(),thrust::plus<rbmd::Real>());

   reduce_by_key(_psample_key.begin(), _psample_key.end(),
    rhok_image_atom.begin(),psamplekey_out.begin(), rhok_image_redue.begin(),
    thrust::equal_to<rbmd::Id>(),thrust::plus<rbmd::Real>());

  //energy
if ("yes" == _energy_rbe_flag ){
  //charge self energy//
  ComputeSelfEnergy(alpha,qqr2e,_e_self_energy);

  //kspace energy
  ComputeKspaceEnergy(box, num_atoms, kmax_array,
      alpha, qqr2e ,_e_kspace);
  _e_kspace = _e_kspace +_e_self_energy;
}

}

void LJCutCoulKspace::ComputeRBE()
{
  auto start = std::chrono::high_resolution_clock::now();
  //
  auto num_atoms = *(_structure_info_data->_num_atoms);
  _rhok_real_redue.resize(_RBE_P);
  _rhok_image_redue.resize(_RBE_P);
  ComputeChargeStructureFactorRBE(*_box, num_atoms, _kmax_array,
      _alpha,_RBE_P,_qqr2e,_rhok_real_redue,_rhok_image_redue);

   //RBE Force
  op::ComputeRBEForceOp<device::DEVICE_GPU>()(
        *_box,num_atoms, _RBE_P,_alpha,_qqr2e,
        thrust::raw_pointer_cast(_rhok_real_redue.data()),
        thrust::raw_pointer_cast(_rhok_image_redue.data()),
        thrust::raw_pointer_cast(_device_data->_d_charge.data()),
        thrust::raw_pointer_cast(_P_Sample_x.data()),
        raw_pointer_cast(_P_Sample_y.data()),
        raw_pointer_cast(_P_Sample_z.data()),
        thrust::raw_pointer_cast(_device_data->_d_px.data()),
        thrust::raw_pointer_cast(_device_data->_d_py.data()),
        thrust::raw_pointer_cast(_device_data->_d_pz.data()),
        thrust::raw_pointer_cast(_device_data->_d_force_kspace_x.data()),
        thrust::raw_pointer_cast(_device_data->_d_force_kspace_y.data()),
        thrust::raw_pointer_cast(_device_data->_d_force_kspace_z.data()),
        thrust::raw_pointer_cast(_device_data->_d_flat_virial_kspace.data()));

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration = end - start;
  TimingStatistics::Instance().record("Long-Range",duration.count());

  //sum virial_kspace on host
  ReduceVirial(num_atoms,_device_data->_d_flat_virial_kspace,
_device_data->_d_virial_kspace);
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

void LJCutCoulKspace::ComputeSelfEnergy(
  rbmd::Real  alpha,
  rbmd::Real  qqr2e,
  rbmd::Real& ave_self_energy)
{
  //compute self_energy
  auto num_atoms = *(_structure_info_data->_num_atoms);
  thrust::device_vector<rbmd::Real> sq_charge;
  sq_charge.resize(num_atoms);

  op::SqchargeOp<device::DEVICE_GPU>()(num_atoms,
    thrust::raw_pointer_cast(_device_data->_d_charge.data()),
    thrust::raw_pointer_cast(sq_charge.data()));

  rbmd::Real sum_sq_charge = thrust::reduce(sq_charge.begin(),
    sq_charge.end(), 0.0f, thrust::plus<rbmd::Real>());
  rbmd::Real total_self_energy = qqr2e * (- sqrt(alpha / M_PI) *sum_sq_charge);

  ave_self_energy =  total_self_energy / num_atoms;
}

void LJCutCoulKspace::ComputeKspaceEnergy(
    Box box,
    rbmd::Id num_atoms,
    Int3 kmax_array,
    rbmd::Real alpha,
    rbmd::Real qqr2e,
    rbmd::Real&  ave_ekspace)
{
    thrust::device_vector<rbmd::Real> density_real_atom;
    thrust::device_vector<rbmd::Real> density_imag_atom;
    density_real_atom.resize(num_atoms);
    density_imag_atom.resize(num_atoms);

    rbmd::Real total_energy_ewald = 0;
    for (rbmd::Id i = -INT_DATA(kmax_array)[0]; i <= INT_DATA(kmax_array)[0]; i++)
    {
        for (rbmd::Id j = -INT_DATA(kmax_array)[1]; j <= INT_DATA(kmax_array)[1]; j++)
        {
            for (rbmd::Id k = -INT_DATA(kmax_array)[2]; k <= INT_DATA(kmax_array)[2]; k++)
            {
                if (!(i == 0 && j == 0 && k == 0))
                {
                    Real3 K = make_Real3(rbmd::Real(2 * M_PI * i / box._length[0]),
                                           rbmd::Real(2 * M_PI * j / box._length[1]),
                                           rbmd::Real(2 * M_PI * k / box._length[2]));
                    rbmd::Real Range_K = SQRT(K.x * K.x + K.y * K.y + K.z * K.z);
                    rbmd::Real Range_K2 = Range_K*Range_K;
                    rbmd::Real alpha_inv =  1 / alpha;

                    op::ComputeChargeStructureFactorOp<device::DEVICE_GPU>()(
                        num_atoms, K,
                        thrust::raw_pointer_cast(_device_data->_d_charge.data()),
                        thrust::raw_pointer_cast(_device_data->_d_px.data()),
                        thrust::raw_pointer_cast(_device_data->_d_py.data()),
                        thrust::raw_pointer_cast(_device_data->_d_pz.data()),
                        thrust::raw_pointer_cast(density_real_atom.data()),
                        thrust::raw_pointer_cast(density_imag_atom.data()));

                    rbmd::Real value_Re = thrust::reduce(density_real_atom.begin(), density_real_atom.end(), 0.0f, thrust::plus<rbmd::Real>());
                    rbmd::Real value_Im = thrust::reduce(density_imag_atom.begin(), density_imag_atom.end(), 0.0f, thrust::plus<rbmd::Real>());
                    rbmd::Real Range_density2 = POW(value_Re, 2.0) + POW(value_Im, 2.0);

                    total_energy_ewald +=
                      EXP(-0.25 * Range_K2 * alpha_inv) * Range_density2 / Range_K2;
                }
            }
        }
    }

  rbmd::Real volume = box._length[0] * box._length[1]*box._length[2];
  total_energy_ewald = qqr2e * (2 * M_PI / volume) * total_energy_ewald;
  ave_ekspace = total_energy_ewald / num_atoms;
}

void LJCutCoulKspace::EvaluatePotentialenergy()
{
  _e_pe_rbl = _e_vdwl_rbl + _e_coul_rbl +_e_kspace;
  //test_ave_pe_rbl = _ave_pe_rbl;

  _e_pe = _e_vdwl+ _e_coul +_e_kspace;
  //test_ave_pe = _ave_pe;

  ThermoStats::Instance().AddThermoData("total-potential-energy",_e_pe);
  
  //out
  std::ofstream outfile("thermo.txt", std::ios::app);
  if (outfile.tellp() == 0) {
    outfile << "step e_vdwl e_coul e_kspace e_pe" << std::endl;
  }
  outfile << test_current_step << " " << _e_vdwl  << " "<< _e_coul <<" "
  << _e_kspace  << " " << _e_pe<< std::endl;
  outfile.close();
}

void LJCutCoulKspace::ComputeQsqSum()
{
  auto num_atoms = *(_structure_info_data->_num_atoms);
  thrust::device_vector<rbmd::Real> sq_charge;
  sq_charge.resize(num_atoms);
  op::SqchargeOp<device::DEVICE_GPU>()(num_atoms,
    thrust::raw_pointer_cast(_device_data->_d_charge.data()),
    thrust::raw_pointer_cast(sq_charge.data()));

   _sum_sq_charge = thrust::reduce(sq_charge.begin(),
    sq_charge.end(), 0.0f, thrust::plus<rbmd::Real>());
  _sum_sq_charge = _qqr2e* _sum_sq_charge;
}

rbmd::Real LJCutCoulKspace::ComputeRMS(rbmd::Id kmax,rbmd::Real box_length,rbmd::Real q2)
{
  auto num_atoms = *(_structure_info_data->_num_atoms);
  rbmd::Real rms = 2.0*q2*_g_ewald/box_length *
    SQRT(1.0/(M_PI*kmax*num_atoms)) *
    EXP(-M_PI*M_PI*kmax*kmax/(_g_ewald*_g_ewald*box_length*box_length));

  return rms;
}

void LJCutCoulKspace::coeffs()
{
  kxvecs.resize(_Kmax3D,0);
  kyvecs.resize(_Kmax3D,0);
  kzvecs.resize(_Kmax3D,0);
  kmax_vec3D.resize(_Kmax3D);
  ug.resize(_Kmax3D,0);
  // std::vector<rbmd::Real> eg_flat(_Kmax3D * 3, 0.0);
  // std::vector<rbmd::Real> vg_flat(_Kmax3D * 6, 0.0);

  eg_flat.resize(_Kmax3D * 3, 0.0);
  vg_flat.resize(_Kmax3D * 6, 0.0);

  auto box =  DataManager::getInstance().getMDData()->_box;
  auto volume = CalculateVolume(*box);
  rbmd::Id k,l,m;
  rbmd::Real sqk,vterm;

  rbmd::Real alpha_inv = 1.0 / _alpha;
  rbmd::Real preu = 4.0*M_PI/volume;

  kcount = 0;

  // (k,0,0), (0,l,0), (0,0,m)

  for (m = 1; m <= _Kmax; m++) {
    sqk = (m*REAL_DATA(_unitk)[0]) * (m*REAL_DATA(_unitk)[0]);
    if (sqk <= _gsqmx) {
      kxvecs[kcount] = m;
      kyvecs[kcount] = 0;
      kzvecs[kcount] = 0;
      ug[kcount] = preu*EXP(-0.25*sqk*alpha_inv)/sqk;
      eg_flat[kcount * 3 + 0]= 2.0*REAL_DATA(_unitk)[0]*m*ug[kcount];
      eg_flat[kcount * 3 + 1]  = 0.0;
      eg_flat[kcount * 3 + 2]  = 0.0;
      vterm = -2.0*(1.0/sqk + 0.25*alpha_inv);
      vg_flat[kcount * 3 + 0] = 1.0 + vterm*(REAL_DATA(_unitk)[0]*m)*(REAL_DATA(_unitk)[0]*m);
      vg_flat[kcount * 3 + 1] = 1.0;
      vg_flat[kcount * 3 + 2] = 1.0;
      vg_flat[kcount * 3 + 3] = 0.0;
      vg_flat[kcount * 3 + 4] = 0.0;
      vg_flat[kcount * 3 + 5] = 0.0;
      kcount++;
    }
    sqk = (m*REAL_DATA(_unitk)[1]) * (m*REAL_DATA(_unitk)[1]);
    if (sqk <= _gsqmx) {
      kxvecs[kcount] = 0;
      kyvecs[kcount] = m;
      kzvecs[kcount] = 0;
      ug[kcount] = preu*EXP(-0.25*sqk*alpha_inv)/sqk;
      eg_flat[kcount * 3 + 0] = 0.0;
      eg_flat[kcount * 3 + 1] = 2.0*REAL_DATA(_unitk)[1]*m*ug[kcount];
      eg_flat[kcount * 3 + 2] = 0.0;
      vterm = -2.0*(1.0/sqk + 0.25*alpha_inv);
      vg_flat[kcount * 3 + 0] = 1.0;
      vg_flat[kcount * 3 + 1] = 1.0 + vterm*(REAL_DATA(_unitk)[1]*m)*(REAL_DATA(_unitk)[1]*m);
      vg_flat[kcount * 3 + 2] = 1.0;
      vg_flat[kcount * 3 + 3] = 0.0;
      vg_flat[kcount * 3 + 4] = 0.0;
      vg_flat[kcount * 3 + 5] = 0.0;
      kcount++;
    }
    sqk = (m*REAL_DATA(_unitk)[2]) * (m*REAL_DATA(_unitk)[2]);
    if (sqk <= _gsqmx) {
      kxvecs[kcount] = 0;
      kyvecs[kcount] = 0;
      kzvecs[kcount] = m;
      ug[kcount] = preu*EXP(-0.25*sqk*alpha_inv)/sqk;
      eg_flat[kcount * 3 + 0] = 0.0;
      eg_flat[kcount * 3 + 1] = 0.0;
      eg_flat[kcount * 3 + 2] = 2.0*REAL_DATA(_unitk)[2]*m*ug[kcount];
      vterm = -2.0*(1.0/sqk + 0.25*alpha_inv);
      vg_flat[kcount * 3 + 0] = 1.0;
      vg_flat[kcount * 3 + 1] = 1.0;
      vg_flat[kcount * 3 + 2] = 1.0 + vterm*(REAL_DATA(_unitk)[2]*m)*(REAL_DATA(_unitk)[2]*m);
      vg_flat[kcount * 3 + 3] = 0.0;
      vg_flat[kcount * 3 + 4] = 0.0;
      vg_flat[kcount * 3 + 5] = 0.0;
      kcount++;
    }
  }
  std::cout << "kcount11: " << kcount  <<std::endl;
  // 1 = (k,l,0), 2 = (k,-l,0)

  for (k = 1; k <= kmax_x; k++) {
    for (l = 1; l <= kmax_y; l++) {
      sqk = (REAL_DATA(_unitk)[0]*k) * (REAL_DATA(_unitk)[0]*k) + (REAL_DATA(_unitk)[1]*l) * (REAL_DATA(_unitk)[1]*l);
      if (sqk <= _gsqmx) {
        kxvecs[kcount] = k;
        kyvecs[kcount] = l;
        kzvecs[kcount] = 0;
        ug[kcount] = preu*EXP(-0.25*sqk*alpha_inv)/sqk;
        eg_flat[kcount * 3 + 0] = 2.0*REAL_DATA(_unitk)[0]*k*ug[kcount];
        eg_flat[kcount * 3 + 1]= 2.0*REAL_DATA(_unitk)[1]*l*ug[kcount];
        eg_flat[kcount * 3 + 2] = 0.0;
        vterm = -2.0*(1.0/sqk + 0.25*alpha_inv);
        vg_flat[kcount * 3 + 0] = 1.0 + vterm*(REAL_DATA(_unitk)[0]*k)*(REAL_DATA(_unitk)[0]*k);
        vg_flat[kcount * 3 + 1] = 1.0 + vterm*(REAL_DATA(_unitk)[1]*l)*(REAL_DATA(_unitk)[1]*l);
        vg_flat[kcount * 3 + 2] = 1.0;
        vg_flat[kcount * 3 + 3] = vterm*REAL_DATA(_unitk)[0]*k*REAL_DATA(_unitk)[1]*l;
        vg_flat[kcount * 3 + 4] = 0.0;
        vg_flat[kcount * 3 + 5] = 0.0;
        kcount++;

        kxvecs[kcount] = k;
        kyvecs[kcount] = -l;
        kzvecs[kcount] = 0;
        ug[kcount] = preu*EXP(-0.25*sqk*alpha_inv)/sqk;
        eg_flat[kcount * 3 + 0] = 2.0*REAL_DATA(_unitk)[0]*k*ug[kcount];
        eg_flat[kcount * 3 + 1]= -2.0*REAL_DATA(_unitk)[1]*l*ug[kcount];
        eg_flat[kcount * 3 + 2] = 0.0;
        vg_flat[kcount * 3 + 0] = 1.0 + vterm*(REAL_DATA(_unitk)[0]*k)*(REAL_DATA(_unitk)[0]*k);
        vg_flat[kcount * 3 + 1] = 1.0 + vterm*(REAL_DATA(_unitk)[1]*l)*(REAL_DATA(_unitk)[1]*l);
        vg_flat[kcount * 3 + 2] = 1.0;
        vg_flat[kcount * 3 + 3] = -vterm*REAL_DATA(_unitk)[0]*k*REAL_DATA(_unitk)[1]*l;
        vg_flat[kcount * 3 + 4] = 0.0;
        vg_flat[kcount * 3 + 5] = 0.0;
        kcount++;
      }
    }
  }

  // 1 = (0,l,m), 2 = (0,l,-m)

  for (l = 1; l <= kmax_y; l++) {
    for (m = 1; m <= kmax_z; m++) {
      sqk = (REAL_DATA(_unitk)[1]*l) * (REAL_DATA(_unitk)[1]*l) + (REAL_DATA(_unitk)[2]*m) * (REAL_DATA(_unitk)[2]*m);
      if (sqk <= _gsqmx) {
        kxvecs[kcount] = 0;
        kyvecs[kcount] = l;
        kzvecs[kcount] = m;
        ug[kcount] = preu*EXP(-0.25*sqk*alpha_inv)/sqk;
        eg_flat[kcount * 3 + 0] =  0.0;
        eg_flat[kcount * 3 + 1]=  2.0*REAL_DATA(_unitk)[1]*l*ug[kcount];
        eg_flat[kcount * 3 + 2] =  2.0*REAL_DATA(_unitk)[2]*m*ug[kcount];
        vterm = -2.0*(1.0/sqk + 0.25*alpha_inv);
        vg_flat[kcount * 3 + 0] = 1.0;
        vg_flat[kcount * 3 + 1] = 1.0 + vterm*(REAL_DATA(_unitk)[1]*l)*(REAL_DATA(_unitk)[1]*l);
        vg_flat[kcount * 3 + 2] = 1.0 + vterm*(REAL_DATA(_unitk)[2]*m)*(REAL_DATA(_unitk)[2]*m);
        vg_flat[kcount * 3 + 3] = 0.0;
        vg_flat[kcount * 3 + 4] = 0.0;
        vg_flat[kcount * 3 + 5] = vterm*REAL_DATA(_unitk)[1]*l*REAL_DATA(_unitk)[2]*m;
        kcount++;

        kxvecs[kcount] = 0;
        kyvecs[kcount] = l;
        kzvecs[kcount] = -m;
        ug[kcount] = preu*EXP(-0.25*sqk*alpha_inv)/sqk;
        eg_flat[kcount * 3 + 0] =  0.0;
        eg_flat[kcount * 3 + 1]=  2.0*REAL_DATA(_unitk)[1]*l*ug[kcount];
        eg_flat[kcount * 3 + 2] = -2.0*REAL_DATA(_unitk)[2]*m*ug[kcount];
        vg_flat[kcount * 3 + 0] = 1.0;
        vg_flat[kcount * 3 + 1] = 1.0 + vterm*(REAL_DATA(_unitk)[1]*l)*(REAL_DATA(_unitk)[1]*l);
        vg_flat[kcount * 3 + 2] = 1.0 + vterm*(REAL_DATA(_unitk)[2]*m)*(REAL_DATA(_unitk)[2]*m);
        vg_flat[kcount * 3 + 3] = 0.0;
        vg_flat[kcount * 3 + 4] = 0.0;
        vg_flat[kcount * 3 + 5] = -vterm*REAL_DATA(_unitk)[1]*l*REAL_DATA(_unitk)[2]*m;
        kcount++;
      }
    }
  }

  // 1 = (k,0,m), 2 = (k,0,-m)

  for (k = 1; k <= kmax_x; k++) {
    for (m = 1; m <= kmax_z; m++) {
      sqk = (REAL_DATA(_unitk)[0]*k) * (REAL_DATA(_unitk)[0]*k) + (REAL_DATA(_unitk)[2]*m) * (REAL_DATA(_unitk)[2]*m);
      if (sqk <= _gsqmx) {
        kxvecs[kcount] = k;
        kyvecs[kcount] = 0;
        kzvecs[kcount] = m;
        ug[kcount] = preu*EXP(-0.25*sqk*alpha_inv)/sqk;
        eg_flat[kcount * 3 + 0] =  2.0*REAL_DATA(_unitk)[0]*k*ug[kcount];
        eg_flat[kcount * 3 + 1]=  0.0;
        eg_flat[kcount * 3 + 2] =  2.0*REAL_DATA(_unitk)[2]*m*ug[kcount];
        vterm = -2.0*(1.0/sqk + 0.25*alpha_inv);
        vg_flat[kcount * 3 + 0] = 1.0 + vterm*(REAL_DATA(_unitk)[0]*k)*(REAL_DATA(_unitk)[0]*k);
        vg_flat[kcount * 3 + 1] = 1.0;
        vg_flat[kcount * 3 + 2] = 1.0 + vterm*(REAL_DATA(_unitk)[2]*m)*(REAL_DATA(_unitk)[2]*m);
        vg_flat[kcount * 3 + 3] = 0.0;
        vg_flat[kcount * 3 + 4] = vterm*REAL_DATA(_unitk)[0]*k*REAL_DATA(_unitk)[2]*m;
        vg_flat[kcount * 3 + 5] = 0.0;
        kcount++;

        kxvecs[kcount] = k;
        kyvecs[kcount] = 0;
        kzvecs[kcount] = -m;
        ug[kcount] = preu*EXP(-0.25*sqk*alpha_inv)/sqk;
        eg_flat[kcount * 3 + 0] =  2.0*REAL_DATA(_unitk)[0]*k*ug[kcount];
        eg_flat[kcount * 3 + 1]=  0.0;
        eg_flat[kcount * 3 + 2] = -2.0*REAL_DATA(_unitk)[2]*m*ug[kcount];
        vg_flat[kcount * 3 + 0] = 1.0 + vterm*(REAL_DATA(_unitk)[0]*k)*(REAL_DATA(_unitk)[0]*k);
        vg_flat[kcount * 3 + 1] = 1.0;
        vg_flat[kcount * 3 + 2] = 1.0 + vterm*(REAL_DATA(_unitk)[2]*m)*(REAL_DATA(_unitk)[2]*m);
        vg_flat[kcount * 3 + 3] = 0.0;
        vg_flat[kcount * 3 + 4] = -vterm*REAL_DATA(_unitk)[0]*k*REAL_DATA(_unitk)[2]*m;
        vg_flat[kcount * 3 + 5] = 0.0;
        kcount++;
      }
    }
  }
  std::cout << "kcount22: " << kcount  <<std::endl;
  // 1 = (k,l,m), 2 = (k,-l,m), 3 = (k,l,-m), 4 = (k,-l,-m)

  for (k = 1; k <= kmax_x; k++) {
    for (l = 1; l <= kmax_y; l++) {
      for (m = 1; m <= kmax_z; m++) {
        sqk = (REAL_DATA(_unitk)[0]*k) * (REAL_DATA(_unitk)[0]*k) + (REAL_DATA(_unitk)[1]*l) * (REAL_DATA(_unitk)[1]*l) +
          (REAL_DATA(_unitk)[2]*m) * (REAL_DATA(_unitk)[2]*m);
        if (sqk <= _gsqmx) {
          kxvecs[kcount] = k;
          kyvecs[kcount] = l;
          kzvecs[kcount] = m;
          ug[kcount] = preu*EXP(-0.25*sqk*alpha_inv)/sqk;
          eg_flat[kcount * 3 + 0] = 2.0*REAL_DATA(_unitk)[0]*k*ug[kcount];
          eg_flat[kcount * 3 + 1]= 2.0*REAL_DATA(_unitk)[1]*l*ug[kcount];
          eg_flat[kcount * 3 + 2] = 2.0*REAL_DATA(_unitk)[2]*m*ug[kcount];
          vterm = -2.0*(1.0/sqk + 0.25*alpha_inv);
          vg_flat[kcount * 3 + 0] = 1.0 + vterm*(REAL_DATA(_unitk)[0]*k)*(REAL_DATA(_unitk)[0]*k);
          vg_flat[kcount * 3 + 1] = 1.0 + vterm*(REAL_DATA(_unitk)[1]*l)*(REAL_DATA(_unitk)[1]*l);
          vg_flat[kcount * 3 + 2] = 1.0 + vterm*(REAL_DATA(_unitk)[2]*m)*(REAL_DATA(_unitk)[2]*m);
          vg_flat[kcount * 3 + 3] = vterm*REAL_DATA(_unitk)[0]*k*REAL_DATA(_unitk)[1]*l;
          vg_flat[kcount * 3 + 4] = vterm*REAL_DATA(_unitk)[0]*k*REAL_DATA(_unitk)[2]*m;
          vg_flat[kcount * 3 + 5] = vterm*REAL_DATA(_unitk)[1]*l*REAL_DATA(_unitk)[2]*m;
          kcount++;

          kxvecs[kcount] = k;
          kyvecs[kcount] = -l;
          kzvecs[kcount] = m;
          ug[kcount] = preu*EXP(-0.25*sqk*alpha_inv)/sqk;
          eg_flat[kcount * 3 + 0] = 2.0*REAL_DATA(_unitk)[0]*k*ug[kcount];
          eg_flat[kcount * 3 + 1]= -2.0*REAL_DATA(_unitk)[1]*l*ug[kcount];
          eg_flat[kcount * 3 + 2] = 2.0*REAL_DATA(_unitk)[2]*m*ug[kcount];
          vg_flat[kcount * 3 + 0] = 1.0 + vterm*(REAL_DATA(_unitk)[0]*k)*(REAL_DATA(_unitk)[0]*k);
          vg_flat[kcount * 3 + 1] = 1.0 + vterm*(REAL_DATA(_unitk)[1]*l)*(REAL_DATA(_unitk)[1]*l);
          vg_flat[kcount * 3 + 2] = 1.0 + vterm*(REAL_DATA(_unitk)[2]*m)*(REAL_DATA(_unitk)[2]*m);
          vg_flat[kcount * 3 + 3] = -vterm*REAL_DATA(_unitk)[0]*k*REAL_DATA(_unitk)[1]*l;
          vg_flat[kcount * 3 + 4] = vterm*REAL_DATA(_unitk)[0]*k*REAL_DATA(_unitk)[2]*m;
          vg_flat[kcount * 3 + 5] = -vterm*REAL_DATA(_unitk)[1]*l*REAL_DATA(_unitk)[2]*m;
          kcount++;

          kxvecs[kcount] = k;
          kyvecs[kcount] = l;
          kzvecs[kcount] = -m;
          ug[kcount] = preu*EXP(-0.25*sqk*alpha_inv)/sqk;
          eg_flat[kcount * 3 + 0] = 2.0*REAL_DATA(_unitk)[0]*k*ug[kcount];
          eg_flat[kcount * 3 + 1]= 2.0*REAL_DATA(_unitk)[1]*l*ug[kcount];
          eg_flat[kcount * 3 + 2] = -2.0*REAL_DATA(_unitk)[2]*m*ug[kcount];
          vg_flat[kcount * 3 + 0] = 1.0 + vterm*(REAL_DATA(_unitk)[0]*k)*(REAL_DATA(_unitk)[0]*k);
          vg_flat[kcount * 3 + 1] = 1.0 + vterm*(REAL_DATA(_unitk)[1]*l)*(REAL_DATA(_unitk)[1]*l);
          vg_flat[kcount * 3 + 2] = 1.0 + vterm*(REAL_DATA(_unitk)[2]*m)*(REAL_DATA(_unitk)[2]*m);
          vg_flat[kcount * 3 + 3] = vterm*REAL_DATA(_unitk)[0]*k*REAL_DATA(_unitk)[1]*l;
          vg_flat[kcount * 3 + 4] = -vterm*REAL_DATA(_unitk)[0]*k*REAL_DATA(_unitk)[2]*m;
          vg_flat[kcount * 3 + 5] = -vterm*REAL_DATA(_unitk)[1]*l*REAL_DATA(_unitk)[2]*m;
          kcount++;

          kxvecs[kcount] = k;
          kyvecs[kcount] = -l;
          kzvecs[kcount] = -m;
          ug[kcount] = preu*EXP(-0.25*sqk*alpha_inv)/sqk;
          eg_flat[kcount * 3 + 0] = 2.0*REAL_DATA(_unitk)[0]*k*ug[kcount];
          eg_flat[kcount * 3 + 1]= -2.0*REAL_DATA(_unitk)[1]*l*ug[kcount];
          eg_flat[kcount * 3 + 2] = -2.0*REAL_DATA(_unitk)[2]*m*ug[kcount];
          vg_flat[kcount * 3 + 0] = 1.0 + vterm*(REAL_DATA(_unitk)[0]*k)*(REAL_DATA(_unitk)[0]*k);
          vg_flat[kcount * 3 + 1] = 1.0 + vterm*(REAL_DATA(_unitk)[1]*l)*(REAL_DATA(_unitk)[1]*l);
          vg_flat[kcount * 3 + 2] = 1.0 + vterm*(REAL_DATA(_unitk)[2]*m)*(REAL_DATA(_unitk)[2]*m);
          vg_flat[kcount * 3 + 3] = -vterm*REAL_DATA(_unitk)[0]*k*REAL_DATA(_unitk)[1]*l;
          vg_flat[kcount * 3 + 4] = -vterm*REAL_DATA(_unitk)[0]*k*REAL_DATA(_unitk)[2]*m;
          vg_flat[kcount * 3 + 5] = vterm*REAL_DATA(_unitk)[1]*l*REAL_DATA(_unitk)[2]*m;
          kcount++;
        }
      }
    }
  }
  std::cout << "kcount33: " << kcount  <<std::endl;
  // thrust::copy(eg_flat.begin(),eg_flat.end(),_d_eg_flat.begin());
  // thrust::copy(vg_flat.begin(),vg_flat.end(),_d_vg_flat.begin());


  _d_kxvecs.resize(_Kmax3D);
  _d_kyvecs.resize(_Kmax3D);
  _d_kzvecs.resize(_Kmax3D);
  MEMCPY(thrust::raw_pointer_cast(_d_kxvecs.data()),kxvecs.data(),
    _Kmax3D * sizeof(rbmd::Id),H2D);
  MEMCPY(thrust::raw_pointer_cast(_d_kyvecs.data()),kyvecs.data(),
  _Kmax3D * sizeof(rbmd::Id),H2D);
  MEMCPY(thrust::raw_pointer_cast(_d_kzvecs.data()),kzvecs.data(),
  _Kmax3D * sizeof(rbmd::Id),H2D);

  _d_eg_flat.resize(_Kmax3D * 3, 0.0);
  _d_vg_flat.resize(_Kmax3D * 6, 0.0);
  MEMCPY(thrust::raw_pointer_cast(_d_eg_flat.data()),eg_flat.data(),
    _Kmax3D * 3* sizeof(rbmd::Real),H2D);
  MEMCPY(thrust::raw_pointer_cast(_d_vg_flat.data()),vg_flat.data(),
    _Kmax3D * 3* sizeof(rbmd::Real),H2D);
}

void LJCutCoulKspace::SetKspacePara()
{
  auto box =  DataManager::getInstance().getMDData()->_box;

  REAL_DATA(_unitk)[0] = 2.0*M_PI/box->_length[0];
  REAL_DATA(_unitk)[1] = 2.0*M_PI/box->_length[1];
  REAL_DATA(_unitk)[2] = 2.0*M_PI/box->_length[2];
  auto num_atoms = *(_structure_info_data->_num_atoms);

  rbmd::Id kmax_x = 1;
  rbmd::Id kmax_y = 1;
  rbmd::Id kmax_z = 1;
  rbmd::Real rms;

  rms = ComputeRMS(kmax_x,box->_length[0],_sum_sq_charge);
  while (rms > _accuracy) {
    kmax_x++;
    rms = ComputeRMS(kmax_x,box->_length[0],_sum_sq_charge);
  }

  rms = ComputeRMS(kmax_y,box->_length[1],_sum_sq_charge);
  while (rms > _accuracy) {
    kmax_y++;
    rms = ComputeRMS(kmax_y,box->_length[1],_sum_sq_charge);
  }

  rms = ComputeRMS(kmax_z,box->_length[2],_sum_sq_charge);
  while (rms > _accuracy) {
    kmax_z++;
    rms = ComputeRMS(kmax_z,box->_length[2],_sum_sq_charge);
  }

  _Kmax = MAX(kmax_x,kmax_y);
  _Kmax = MAX(_Kmax,kmax_z);
  _Kmax3D = 4*_Kmax*_Kmax*_Kmax + 6*_Kmax*_Kmax + 3*_Kmax;
  _kmax_array = {kmax_x,kmax_y,kmax_z};

  std::cout<< "Kmax: "<< _Kmax << ", Kmax3D: "<<  _Kmax3D <<std::endl;

  rbmd::Real gsqxmx = REAL_DATA(_unitk)[0] *REAL_DATA(_unitk)[0] *kmax_x*kmax_x;
  rbmd::Real gsqymx = REAL_DATA(_unitk)[1] *REAL_DATA(_unitk)[1] *kmax_y*kmax_y;
  rbmd::Real gsqzmx = REAL_DATA(_unitk)[2] *REAL_DATA(_unitk)[2] *kmax_z*kmax_z;
  _gsqmx = MAX(gsqxmx,gsqymx);
  _gsqmx = MAX(_gsqmx,gsqzmx);


  auto kmax_read_flag = 0 ;
  rbmd::Id kmax_x_read,kmax_y_read,kmax_z_read;
  if(kmax_read_flag)
  {
    kmax_x = kmax_x_read;
    kmax_y = kmax_y_read;
    kmax_z = kmax_z_read;

    _Kmax = MAX(kmax_x,kmax_y);
    _Kmax = MAX(_Kmax,kmax_z);
    _Kmax3D = 4*_Kmax*_Kmax*_Kmax + 6*_Kmax*_Kmax + 3*_Kmax;
    _kmax_array = {kmax_x,kmax_y,kmax_z};

    rbmd::Real gsqxmx = REAL_DATA(_unitk)[0] *REAL_DATA(_unitk)[0] *kmax_x*kmax_x;
    rbmd::Real gsqymx = REAL_DATA(_unitk)[1] *REAL_DATA(_unitk)[1] *kmax_y*kmax_y;
    rbmd::Real gsqzmx = REAL_DATA(_unitk)[2] *REAL_DATA(_unitk)[2] *kmax_z*kmax_z;
    _gsqmx = MAX(gsqxmx,gsqymx);
    _gsqmx = MAX(_gsqmx,gsqzmx);
  }
  _gsqmx *= 1.00001;
  //coeffs();//
  std::cout << "gsqmx: " << _gsqmx << " " << "kcount: " << kcount  <<std::endl;
}


void LJCutCoulKspace::ComputeQsf()
{
  auto num_atoms = *(_structure_info_data->_num_atoms);
  _d_cs.resize(num_atoms * 3 *  (2 * _Kmax + 1));
  _d_sn.resize(num_atoms * 3 *  (2 * _Kmax + 1));
  _d_qfactor_real.resize(_Kmax3D);
  _d_qfactor_image.resize(_Kmax3D);

  thrust::fill(_d_qfactor_real.begin(),
  _d_qfactor_real.end(), 0.0f);
  thrust::fill(_d_qfactor_image.begin(),
    _d_qfactor_image.end(), 0.0f);

  op::EikOp<device::DEVICE_GPU>()(
    num_atoms,_gsqmx,_unitk,_Kmax,_kmax_array,
    thrust::raw_pointer_cast(_device_data->_d_px.data()),
 thrust::raw_pointer_cast(_device_data->_d_py.data()),
 thrust::raw_pointer_cast(_device_data->_d_pz.data()),
 thrust::raw_pointer_cast(_device_data->_d_charge.data()),
 thrust::raw_pointer_cast(_d_cs.data()),
 thrust::raw_pointer_cast(_d_sn.data()),
 thrust::raw_pointer_cast(_d_qfactor_real.data()),
 thrust::raw_pointer_cast(_d_qfactor_image.data()));

  std::vector<rbmd::Real> qfactor_real(_Kmax3D);
  std::vector<rbmd::Real> qfactor_imag(_Kmax3D);


  thrust::copy(_d_qfactor_real.begin(),
    _d_qfactor_real.end(), qfactor_real.begin());

  thrust::copy(_d_qfactor_image.begin(),
  _d_qfactor_image.end(), qfactor_imag.begin());


  std::ofstream output_file("output_qfactor.txt");
  for (size_t i = 0; i < qfactor_real.size(); ++i)
  {
    output_file << "i:" << i << " "
    << qfactor_real[i] << " " << qfactor_imag[i]  << std::endl;

  }
}

void LJCutCoulKspace::ComputeQsf_fix()
{
  auto start = std::chrono::high_resolution_clock::now();
  _d_qfactor_real.resize(_Kmax3D);
  _d_qfactor_image.resize(_Kmax3D);

  auto num_atoms = *(_structure_info_data->_num_atoms);
  auto atom_id_to_idx =
  LinkedCellLocator::GetInstance().GetLinkedCell()->_atom_id_to_idx;
  op::EikFixOP<device::DEVICE_GPU>()(
    num_atoms,kcount,_gsqmx,_unitk,
    thrust::raw_pointer_cast(_d_kxvecs.data()),
    thrust::raw_pointer_cast(_d_kyvecs.data()),
    thrust::raw_pointer_cast(_d_kzvecs.data()),
    thrust::raw_pointer_cast(_device_data->_d_atoms_id.data()),
    thrust::raw_pointer_cast(atom_id_to_idx.data()),
    thrust::raw_pointer_cast(_device_data->_d_px.data()),
    thrust::raw_pointer_cast(_device_data->_d_py.data()),
    thrust::raw_pointer_cast(_device_data->_d_pz.data()),
    thrust::raw_pointer_cast(_device_data->_d_charge.data()),
    thrust::raw_pointer_cast(_d_qfactor_real.data()),
    thrust::raw_pointer_cast(_d_qfactor_image.data()));


  auto end = std::chrono::high_resolution_clock::now();

  std::chrono::duration<rbmd::Real> duration = end - start;
  std::cout << "构建qfactor耗时" << duration.count() << "秒" << std::endl;

  std::vector<rbmd::Real> qfactor_real(_Kmax3D);
  std::vector<rbmd::Real> qfactor_image(_Kmax3D);

  thrust::copy(_d_qfactor_real.begin(),
    _d_qfactor_real.end(), qfactor_real.begin());
  thrust::copy(_d_qfactor_image.begin(),
_d_qfactor_image.end(), qfactor_image.begin());


  std::ofstream output_file("output_qfactor.txt");
  for (size_t i = 0; i < qfactor_real.size(); ++i)
  {
    output_file << "i:" << i << " "
    << qfactor_real[i] << " " << qfactor_image[i]  << std::endl;
  }
  output_file.close();


}

void LJCutCoulKspace::ComputeWaveVectors()
{
 // 首先遍历所有波矢组合，计算满足条件的波矢总数
    int total_wavevectors = 0;

    // 1. 处理单分量波矢 (k,0,0), (0,l,0), (0,0,m)
    for(int ic = 0; ic < 3; ic++) {
        int max_idx;
        if(ic == 0) max_idx = kmax_x;
        else if(ic == 1) max_idx = kmax_y;
        else max_idx = kmax_z;

        for(int m = 1; m <= max_idx; m++) {
            double sqk = (m * REAL_DATA(_unitk)[ic]) * (m * REAL_DATA(_unitk)[ic]);
            if(sqk <= _gsqmx) {
                total_wavevectors += 1; // 每个符合条件的单分量波矢计数一次
            }
        }
    }

    // 2. 处理双分量波矢 (k,l,0), (k,-l,0)
    for(int k = 1; k <= kmax_x; k++) {
        for(int l = 1; l <= kmax_y; l++) {
            double sqk = (k * REAL_DATA(_unitk)[0]) * (k * REAL_DATA(_unitk)[0]) +
                         (l * REAL_DATA(_unitk)[1]) * (l * REAL_DATA(_unitk)[1]);
            if(sqk <= _gsqmx) {
                total_wavevectors += 2; // (k,l,0) 和 (k,-l,0) 各计数一次
            }
        }
    }

    // 3. 处理双分量波矢 (0,l,m), (0,l,-m)
    for(int l = 1; l <= kmax_y; l++) {
        for(int m = 1; m <= kmax_z; m++) {
            double sqk = (l * REAL_DATA(_unitk)[1]) * (l * REAL_DATA(_unitk)[1]) +
                         (m * REAL_DATA(_unitk)[2]) * (m * REAL_DATA(_unitk)[2]);
            if(sqk <= _gsqmx) {
                total_wavevectors += 2; // (0,l,m) 和 (0,l,-m) 各计数一次
            }
        }
    }

    // 4. 处理双分量波矢 (k,0,m), (k,0,-m)
    for(int k = 1; k <= kmax_x; k++) {
        for(int m = 1; m <= kmax_z; m++) {
            double sqk = (k * REAL_DATA(_unitk)[0]) * (k * REAL_DATA(_unitk)[0]) +
                         (m * REAL_DATA(_unitk)[2]) * (m * REAL_DATA(_unitk)[2]);
            if(sqk <= _gsqmx) {
                total_wavevectors += 2; // (k,0,m) 和 (k,0,-m) 各计数一次
            }
        }
    }

    // 5. 处理三分量波矢 (k,l,m), (k,-l,m), (k,l,-m), (k,-l,-m)
    for(int k = 1; k <= kmax_x; k++) {
        for(int l = 1; l <= kmax_y; l++) {
            for(int m = 1; m <= kmax_z; m++) {
                double sqk = (k * REAL_DATA(_unitk)[0]) * (k * REAL_DATA(_unitk)[0]) +
                             (l * REAL_DATA(_unitk)[1]) * (l * REAL_DATA(_unitk)[1]) +
                             (m * REAL_DATA(_unitk)[2]) * (m * REAL_DATA(_unitk)[2]);
                if(sqk <= _gsqmx) {
                    total_wavevectors += 4; // 四种组合各计数一次
                }
            }
        }
    }


  ///////
  wavevec_indices_host.resize(3*total_wavevectors);

  // 第二次遍历，填充 wavevec_indices_host
    int n = 0;

    // 1. 填充单分量波矢 (k,0,0), (0,l,0), (0,0,m)
    for(int ic = 0; ic < 3; ic++)
    {
        int max_idx;
        if(ic == 0) max_idx = kmax_x;
        else if(ic == 1) max_idx = kmax_y;
        else max_idx = kmax_z;

        for(int m = 1; m <= max_idx; m++) {
            double sqk = (m * REAL_DATA(_unitk)[ic]) * (m * REAL_DATA(_unitk)[ic]);
            if(sqk <= _gsqmx) {
                int k_val = 0, l_val = 0, m_val = 0;
                if(ic == 0) { k_val = m; l_val = 0; m_val = 0; }
                else if(ic == 1) { k_val = 0; l_val = m; m_val = 0; }
                else { k_val = 0; l_val = 0; m_val = m; }

                wavevec_indices_host[3*n + 0] = k_val;
                wavevec_indices_host[3*n + 1] = l_val;
                wavevec_indices_host[3*n + 2] = m_val;
                n++;
            }
        }
    }

    // 2. 填充双分量波矢 (k,l,0), (k,-l,0)
    for(int k = 1; k <= kmax_x; k++) {
        for(int l = 1; l <= kmax_y; l++) {
            double sqk = (k * REAL_DATA(_unitk)[0]) * (k * REAL_DATA(_unitk)[0]) +
                         (l * REAL_DATA(_unitk)[1]) * (l * REAL_DATA(_unitk)[1]);
            if(sqk <= _gsqmx) {
                // (k, l, 0)
                wavevec_indices_host[3*n + 0] = k;
                wavevec_indices_host[3*n + 1] = l;
                wavevec_indices_host[3*n + 2] = 0;
                n++;

                // (k, -l, 0)
                wavevec_indices_host[3*n + 0] = k;
                wavevec_indices_host[3*n + 1] = -l;
                wavevec_indices_host[3*n + 2] = 0;
                n++;
            }
        }
    }

    // 3. 填充双分量波矢 (0,l,m), (0,l,-m)
    for(int l = 1; l <= kmax_y; l++) {
        for(int m = 1; m <= kmax_z; m++) {
            double sqk = (l * REAL_DATA(_unitk)[1]) * (l * REAL_DATA(_unitk)[1]) +
                         (m * REAL_DATA(_unitk)[2]) * (m * REAL_DATA(_unitk)[2]);
            if(sqk <= _gsqmx) {
                // (0, l, m)
                wavevec_indices_host[3*n + 0] = 0;
                wavevec_indices_host[3*n + 1] = l;
                wavevec_indices_host[3*n + 2] = m;
                n++;

                // (0, l, -m)
                wavevec_indices_host[3*n + 0] = 0;
                wavevec_indices_host[3*n + 1] = l;
                wavevec_indices_host[3*n + 2] = -m;
                n++;
            }
        }
    }

    // 4. 填充双分量波矢 (k,0,m), (k,0,-m)
    for(int k = 1; k <= kmax_x; k++) {
        for(int m = 1; m <= kmax_z; m++) {
            double sqk = (k * REAL_DATA(_unitk)[0]) * (k * REAL_DATA(_unitk)[0]) +
                         (m * REAL_DATA(_unitk)[2]) * (m * REAL_DATA(_unitk)[2]);
            if(sqk <= _gsqmx) {
                // (k, 0, m)
                wavevec_indices_host[3*n + 0] = k;
                wavevec_indices_host[3*n + 1] = 0;
                wavevec_indices_host[3*n + 2] = m;
                n++;

                // (k, 0, -m)
                wavevec_indices_host[3*n + 0] = k;
                wavevec_indices_host[3*n + 1] = 0;
                wavevec_indices_host[3*n + 2] = -m;
                n++;
            }
        }
    }

    // 5. 填充三分量波矢 (k,l,m), (k,-l,m), (k,l,-m), (k,-l,-m)
    for(int k = 1; k <= kmax_x; k++) {
        for(int l = 1; l <= kmax_y; l++) {
            for(int m = 1; m <= kmax_z; m++) {
                double sqk = (k * REAL_DATA(_unitk)[0]) * (k * REAL_DATA(_unitk)[0]) +
                             (l * REAL_DATA(_unitk)[1]) * (l * REAL_DATA(_unitk)[1]) +
                             (m * REAL_DATA(_unitk)[2]) * (m * REAL_DATA(_unitk)[2]);
                if(sqk <= _gsqmx) {
                    // (k, l, m)
                    wavevec_indices_host[3*n + 0] = k;
                    wavevec_indices_host[3*n + 1] = l;
                    wavevec_indices_host[3*n + 2] = m;
                    n++;

                    // (k, -l, m)
                    wavevec_indices_host[3*n + 0] = k;
                    wavevec_indices_host[3*n + 1] = -l;
                    wavevec_indices_host[3*n + 2] = m;
                    n++;

                    // (k, l, -m)
                    wavevec_indices_host[3*n + 0] = k;
                    wavevec_indices_host[3*n + 1] = l;
                    wavevec_indices_host[3*n + 2] = -m;
                    n++;

                    // (k, -l, -m)
                    wavevec_indices_host[3*n + 0] = k;
                    wavevec_indices_host[3*n + 1] = -l;
                    wavevec_indices_host[3*n + 2] = -m;
                    n++;
                }
            }
        }
    }

  // 确认填充的波矢数量与计算的total_wavevectors一致
  if(n != total_wavevectors) {
    fprintf(stderr, "Mismatch in wavevector count: counted=%d, filled=%d\n",
      total_wavevectors, n);
  }
  printf("total_wavevectors: %i     n:  %i\n",total_wavevectors,n);
}

void LJCutCoulKspace::ComputeEwlad_fix()
{
  //
  //ComputeQsf_fix();

  auto box =  DataManager::getInstance().getMDData()->_box;
  ComputeWaveVectors_2(*box,_Kmax);
  //ComputeQsf();

  // auto num_atoms = *(_structure_info_data->_num_atoms);
  // // charge structure factors
  // for (rbmd::Id k_index = 0; k_index < kcount; k_index++)
  // {
  //   rbmd::Id kx = kxvecs[k_index];
  //   rbmd::Id ky = kyvecs[k_index];
  //   rbmd::Id kz = kzvecs[k_index];
  //   Int3 kmax_vec3D ={kx,ky,kz};
  //   op::EwaldForceFixOp<device::DEVICE_GPU>()(
  //     num_atoms,kcount,k_index,_qqr2e,kmax_vec3D,
  //     thrust::raw_pointer_cast(_d_eg_flat.data()),
  //     thrust::raw_pointer_cast(_d_cs.data()),
  //     thrust::raw_pointer_cast(_d_sn.data()),
  //   thrust::raw_pointer_cast(_device_data->_d_charge.data()),
  //   thrust::raw_pointer_cast(_d_qfactor_real.data()),
  //   thrust::raw_pointer_cast(_d_qfactor_image.data()),
  //     thrust::raw_pointer_cast(_device_data->_d_force_kspace_x.data()),
  //     thrust::raw_pointer_cast(_device_data->_d_force_kspace_y.data()),
  //     thrust::raw_pointer_cast(_device_data->_d_force_kspace_z.data()));
  //   }
  //
  //


  // std::cout << "Ewlad_fix: " << kcount  <<std::endl;
  //
  // std::vector<rbmd::Real> h_force_kspace_x(num_atoms);
  // std::vector<rbmd::Real> h_force_kspace_y(num_atoms);
  // std::vector<rbmd::Real> h_force_kspace_z(num_atoms);
  //
  // thrust::copy(_device_data->_d_force_kspace_x.begin(),
  //   _device_data->_d_force_kspace_x.end(), h_force_kspace_x.begin());
  // thrust::copy(_device_data->_d_force_kspace_y.begin(),
  // _device_data->_d_force_kspace_y.end(), h_force_kspace_y.begin());
  // thrust::copy(_device_data->_d_force_kspace_z.begin(),
  // _device_data->_d_force_kspace_z.end(), h_force_kspace_z.begin());
  //
  // std::ofstream output_file("output_force_kspace111.txt");
  // for (size_t i = 0; i < h_force_kspace_x.size(); ++i)
  // {
  //   output_file << "i:" << i << " "
  //   << h_force_kspace_x[i] << " " << h_force_kspace_y[i]  << " " << h_force_kspace_z[i]
  //   << std::endl;
  // }
  // output_file.close();
}

void LJCutCoulKspace::ComputeWaveVectors_2( Box box,rbmd::Id Kmax)
{
  kxvecs_R.resize(_num_k);
  kyvecs_R.resize(_num_k);
  kzvecs_R.resize(_num_k);
  rbmd::Id count_K = 0;
  for (rbmd::Id i = -Kmax; i <= Kmax; i++)
  {
    for (rbmd::Id j = -Kmax; j <= Kmax; j++)
    {
      for (rbmd::Id k = -Kmax; k <= Kmax; k++)
      {
        if (!(i == 0 && j == 0 && k == 0))
        {
          kxvecs_R[count_K]= 2 * M_PI * i /  box._length[0];
          kyvecs_R[count_K]= 2 * M_PI * j /  box._length[1];
          kzvecs_R[count_K]= 2 * M_PI * k /  box._length[2];
          count_K++;
        }
      }
    }
  }

  _d_kxvecs_R.resize(_num_k);
  _d_kyvecs_R.resize(_num_k);
  _d_kzvecs_R.resize(_num_k);
  thrust::copy(kxvecs_R.begin(),kxvecs_R.end(), _d_kxvecs_R.begin());
  thrust::copy(kyvecs_R.begin(),kyvecs_R.end(), _d_kyvecs_R.begin());
  thrust::copy(kzvecs_R.begin(),kzvecs_R.end(), _d_kzvecs_R.begin());


  auto start = std::chrono::high_resolution_clock::now();
  _d_qfactor_real.resize(_num_k);
  _d_qfactor_image.resize(_num_k);

  auto num_atoms = *(_structure_info_data->_num_atoms);
  auto atom_id_to_idx =
  LinkedCellLocator::GetInstance().GetLinkedCell()->_atom_id_to_idx;
  op::EikFix3OP<device::DEVICE_GPU>()(
    num_atoms,_num_k,
    thrust::raw_pointer_cast(_d_kxvecs_R.data()),
    thrust::raw_pointer_cast(_d_kyvecs_R.data()),
    thrust::raw_pointer_cast(_d_kzvecs_R.data()),
    thrust::raw_pointer_cast(_device_data->_d_atoms_id.data()),
    thrust::raw_pointer_cast(atom_id_to_idx.data()),
    thrust::raw_pointer_cast(_device_data->_d_px.data()),
    thrust::raw_pointer_cast(_device_data->_d_py.data()),
    thrust::raw_pointer_cast(_device_data->_d_pz.data()),
    thrust::raw_pointer_cast(_device_data->_d_charge.data()),
    thrust::raw_pointer_cast(_d_qfactor_real.data()),
    thrust::raw_pointer_cast(_d_qfactor_image.data()));
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration = end - start;
  std::cout << "构建优化 qactor耗时" << duration.count() << "秒" << std::endl;

  std::vector<rbmd::Real> qfactor_real(_num_k);
  std::vector<rbmd::Real> qfactor_image(_num_k);

  thrust::copy(_d_qfactor_real.begin(),
    _d_qfactor_real.end(), qfactor_real.begin());
  thrust::copy(_d_qfactor_image.begin(),
_d_qfactor_image.end(), qfactor_image.begin());


  std::ofstream output_file("output_qfactor.txt");
  for (size_t i = 0; i < qfactor_real.size(); ++i)
  {
    output_file << "i:" << i << " "
    << qfactor_real[i] << " " << qfactor_image[i]  << std::endl;
  }
  output_file.close();
}