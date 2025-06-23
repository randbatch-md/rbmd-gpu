#include "cvff.h"

#include <output/include/Logger.hpp>

#include "../../common/device_types.h"
#include "../../common/rbmd_define.h"
#include "../../common/types.h"
#include "../common/RBEPSample.h"
#include "../common/unit_factor.h"
#include "force_op/force_op.h"
#include "cvff_op/cvff_op.h"
#include "lj_cut_coul_kspace_op/lj_cut_coul_kspace_op.h"
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


//
//   //automatically compute kmax
//   SetKspacePara(); //kmax
//   std::cout <<"g_ewald: " << _g_ewald << ", alpha: "<< _alpha
//     << ", Kmax: " << _kmax_array.x << " " << _kmax_array.y << " "
//     << _kmax_array.z << std::endl;

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

  //coulomb
  _coulomb_type = "NULL"; // default
  if (config->PathExists({"hyper_parameters", "coulomb"}))
  {
    //accuracy
    if (config->PathExists({"hyper_parameters", "coulomb" ,"accuracy"})) {
      _accuracy = DataManager::getInstance().getConfigData()->Get<rbmd::Real>(
"accuracy", "hyper_parameters", "coulomb");

      auto box =  DataManager::getInstance().getMDData()->_box;
      auto volue = CalculateVolume(*box);
      auto num_atoms = *(_structure_info_data->_num_atoms);
      //  sum q_sq
      ComputeQsqSum(); //q2

      //compute g_ewald
      _g_ewald = _accuracy*SQRT(num_atoms*_cut_off*volue) / (2.0*_sum_sq_charge);
      if (_g_ewald >= 1.0) _g_ewald = (1.35 - 0.15*LOG(_accuracy))/_cut_off;
      else _g_ewald = SQRT(-LOG(_g_ewald)) / _cut_off;
      _alpha = _g_ewald*_g_ewald;
      std::cout << "accuracy : " <<_accuracy  << ",  alpha= " <<  _alpha <<std::endl;
    }

    //alpha
    if (config->PathExists({"hyper_parameters", "coulomb" ,"alpha"})) {
      _alpha = config->Get<rbmd::Real>("alpha", "hyper_parameters", "coulomb");

      auto accuracy_test = ERFC(_cut_off * SQRT(_alpha));
      //std::cout << "accuracy_test= " <<  accuracy_test <<std::endl;
    }

    //Kmax
    auto Kmax =config->GetArray<rbmd::Id>("kmax", "hyper_parameters", "coulomb");
    _kmax_array.x = Kmax[0];
    _kmax_array.y = Kmax[1];
    _kmax_array.z = Kmax[2];
    _num_k =  (2*_kmax_array.x +1)  * (2*_kmax_array.y +1)
              * (2*_kmax_array.z +1) - 1;

    //RBE
    _coulomb_type = config->Get<std::string>("type", "hyper_parameters", "coulomb");
    if("RBE" == _coulomb_type) {
      _RBE_P = config->Get<rbmd::Id>("coulomb_sample_num", "hyper_parameters", "coulomb");
      GetPsampleKey();

      //
      bool energy_rbe_flag = config->PathExists({"hyper_parameters", "coulomb" ,"energy_rbe_flag"});
      if (energy_rbe_flag) {
        _energy_rbe_flag = config->Get<std::string>("energy_rbe_flag", "hyper_parameters", "coulomb");
      }
      else {
        Logger::Instance().error( "\033[31m When using RBE for the coulomb type, "
                     "the key 'energy_rbe_flag' must be defined.\033[0m");
        exit(EXIT_FAILURE); //
      }
    }
  }
}

void CVFF::Execute() {
  const auto& config = DataManager::getInstance().getConfigData();

  ComputeLJCutCoulForce();
  if(config->PathExists({"hyper_parameters", "coulomb"})) {
    ComputeKspaceForce();
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

void CVFF::ComputeKspaceForce()
{
  if("RBE" == _coulomb_type)
  {
      ComputeRBE();
  }
  else
  {
     ComputeEwlad();
  }
  ThermoStats::Instance().AddThermoData("kspace",_e_kspace);
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

//   auto num_atoms = *(_structure_info_data->_num_atoms);
//   auto atom_id_to_idx =
// LinkedCellLocator::GetInstance().GetLinkedCell()->_atom_id_to_idx;
//
//   thrust::host_vector<rbmd::Real> h_fx(num_atoms);
//   thrust::host_vector<rbmd::Real> h_fy(num_atoms);
//   thrust::host_vector<rbmd::Real> h_fz(num_atoms);
//
//
//   thrust::copy(_device_data->_d_fx.begin(),_device_data->_d_fx.end(),
//     h_fx.begin());
//   thrust::copy(_device_data->_d_fy.begin(),_device_data->_d_fy.end(),
//     h_fy.begin());
//   thrust::copy(_device_data->_d_fz.begin(),_device_data->_d_fz.end(),
//   h_fz.begin());

  // std::ofstream fx("f_cvff.txt");
  // if (fx.is_open()) {
  //   for (rbmd::Id i = 0; i < num_atoms; ++i) {
  //     auto idx = atom_id_to_idx[i];
  //     fx << i  << " " << h_fx[idx] << " " << h_fy[idx]
  //       << " " << h_fz[idx]  << "\n";
  //   }
  //   fx.close();
  // }
}

void CVFF::ComputeChargeStructureFactorEwald(
    Box box,
    rbmd::Id num_atoms,
    Int3 Kmax_array,
    rbmd::Real alpha,
    rbmd::Real qqr2e,
    thrust::host_vector<rbmd::Real> value_Re_array,
    thrust::host_vector<rbmd::Real> value_Im_array)
{
    //thrust::fill(density_real.begin(), density_real.end(), 0.0f);
    //thrust::fill(density_imag.begin(), density_imag.end(), 0.0f);
    thrust::device_vector<rbmd::Real> density_real_atom;
    thrust::device_vector<rbmd::Real> density_imag_atom;
    density_real_atom.resize(num_atoms);
    density_imag_atom.resize(num_atoms);

    rbmd::Real total_energy_kspace= 0;
    rbmd::Id index = 0;
    for (rbmd::Id i = -INT_DATA(Kmax_array)[0]; i <= INT_DATA(Kmax_array)[0]; i++)
    {
        for (rbmd::Id j = -INT_DATA(Kmax_array)[1]; j <= INT_DATA(Kmax_array)[1]; j++)
        {
            for (rbmd::Id k = -INT_DATA(Kmax_array)[2]; k <= INT_DATA(Kmax_array)[2]; k++)
            {
                if (!(i == 0 && j == 0 && k == 0))
                {
                    Real3 K = make_Real3(rbmd::Real(2 * M_PI * i / box._length[0]),
                                           rbmd::Real(2 * M_PI * j / box._length[1]),
                                           rbmd::Real(2 * M_PI * k / box._length[2]));
                    rbmd::Real Range_K = SQRT(K.x * K.x + K.y * K.y + K.z * K.z);
                    rbmd::Real Range_K2 = Range_K*Range_K;
                    rbmd::Real alpha_inv =  1 / alpha;

                    op::ComputeChargeStructureFactorOp<device::DEVICE_GPU>()(num_atoms, K,
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

void CVFF::ComputeEwlad()
{
  auto start = std::chrono::high_resolution_clock::now();

  auto num_atoms = *(_structure_info_data->_num_atoms);

  //compute charge structure factor
  thrust::host_vector<rbmd::Real> h_Re_array(_num_k);
  thrust::host_vector<rbmd::Real> h_Im_array(_num_k);
  ComputeChargeStructureFactorEwald(*_box, num_atoms, _kmax_array,
    _alpha,_qqr2e, h_Re_array, h_Im_array);

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


void CVFF::RBEInit(Box box,rbmd::Real alpha,rbmd::Id RBE_P)
{
  Real3 sigma = { rbmd::Real((SQRT(alpha / 2.0) * box._length[0]/M_PI)),
                  rbmd::Real((SQRT(alpha / 2.0) * box._length[1]/M_PI)),
                  rbmd::Real((SQRT(alpha / 2.0) * box._length[2]/M_PI))};
  auto random = true;
  RBEPSAMPLE rbe_presolve_psample = { alpha, box, RBE_P, random};
  thrust::host_vector<rbmd::Real> h_P_Sample_x(RBE_P);
  thrust::host_vector<rbmd::Real> h_P_Sample_y(RBE_P);
  thrust::host_vector<rbmd::Real> h_P_Sample_z(RBE_P);

  // TODO Reconstruct using random number generator!
  rbe_presolve_psample.Fetch_P_Sample(0.0, sigma,
    thrust::raw_pointer_cast(h_P_Sample_x.data()),
    thrust::raw_pointer_cast(h_P_Sample_y.data()),
    thrust::raw_pointer_cast(h_P_Sample_z.data()));
  _P_Sample_x = h_P_Sample_x;
  _P_Sample_y = h_P_Sample_y;
  _P_Sample_z = h_P_Sample_z;
}

void CVFF::GetPsampleKey()
{
  //psample index key
  auto num_atoms = *(_structure_info_data->_num_atoms);
  _psample_key.resize(num_atoms * _RBE_P);
  op::GenerateIndexArrayOp<device::DEVICE_GPU>()(num_atoms,_RBE_P,
    thrust::raw_pointer_cast(_psample_key.data()));
}

void CVFF::ComputeChargeStructureFactorRBE(
   Box box,
   rbmd::Id num_atoms,
   Int3 kmax_array,
   rbmd::Real alpha,
   rbmd::Id RBE_P,
   rbmd::Real qqr2e,
   thrust::device_vector<rbmd::Real> rhok_real_redue,
   thrust::device_vector<rbmd::Real> rhok_image_redue)
{
  //get P_Sample at each step
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
  if ("yes" == _energy_rbe_flag )
  {
    //charge self energy//
    ComputeSelfEnergy(alpha,qqr2e,_e_self_energy);

    //kspace energy
    ComputeKspaceEnergy(box, num_atoms, kmax_array,
         alpha, qqr2e ,_e_kspace);
    _e_kspace = _e_kspace +_e_self_energy;
  }
}

void CVFF::ComputeRBE()
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

void CVFF::ComputeSelfEnergy(
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

void CVFF::ComputeKspaceEnergy(
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

                    op::ComputeChargeStructureFactorOp<device::DEVICE_GPU>()(num_atoms, K,
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

void CVFF::ComputeQsqSum()
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

rbmd::Real CVFF::ComputeRMS(rbmd::Id kmax,rbmd::Real box_length,rbmd::Real q2)
{
  auto num_atoms = *(_structure_info_data->_num_atoms);
  rbmd::Real rms = 2.0*q2*_g_ewald/box_length *
    SQRT(1.0/(M_PI*kmax*num_atoms)) *
    EXP(-M_PI*M_PI*kmax*kmax/(_g_ewald*_g_ewald*box_length*box_length));

  return rms;
}

void CVFF::SetKspacePara()
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

  //std::cout << "gsqmx: " << _gsqmx << " " << "kcount: " << kcount  <<std::endl;
}


