#include "kspace_calculator.h"

#include <chrono>
#include <cmath>

#include "../common/unit_factor.h"
#include "common/RBEPSample.h"
#include "common/thermo_stats.hpp"
#include "common/timing_statistics.hpp"
#include "data_manager.h"
#include "force.h"
#include "kspace_calculator_op/kspace_calculator_op.h"
#include "neighbor_list/include/linked_cell/linked_cell_locator.h"
#include "output/include/Logger.hpp"

#define SMALL 0.00001
KSpaceCalculator::KSpaceCalculator()
{

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
}

void KSpaceCalculator::Init()
{
    const auto& config = DataManager::getInstance().getConfigData();
    _cut_off = config->Get<rbmd::Real>("cut_off", "hyper_parameters", "neighbor");

  //coulomb
  _coulomb_type = "NULL"; // default
  if (config->PathExists({"hyper_parameters", "coulomb"}))
  {
    ComputeQsqSum(); //q2

    //accuracy
    if (config->PathExists({"hyper_parameters", "coulomb" ,"accuracy"})) {
      auto relative_accuracy_user = DataManager::getInstance().getConfigData()->Get<rbmd::Real>(
"accuracy", "hyper_parameters", "coulomb");
      auto qelectron = 1; //e
      auto angstrom = 1; //A
      double two_charge_force = _qqr2e* (qelectron*qelectron) /(angstrom*angstrom);
      _accuracy = relative_accuracy_user * two_charge_force;

      auto box =  DataManager::getInstance().getMDData()->_box;
      double volume = CalculateVolume(*box);
      auto num_atoms = *(_structure_info_data->_num_atoms);
      //  sum q_sq
      ComputeQsqSum(); //q2

      // 详细调试输出
      std::cout << "=== Ewald Parameter Debug ===" << std::endl;
      std::cout << "relative_accuracy_user: " << relative_accuracy_user << std::endl;
      std::cout << "two_charge_force: " << two_charge_force << std::endl;
      std::cout << "_accuracy: " << _accuracy << std::endl;
      std::cout << "num_atoms: " << num_atoms << std::endl;
      std::cout << "_cut_off: " << _cut_off << std::endl;
      std::cout << "volume: " << volume << std::endl;
      std::cout << "_sum_sq_charge: " << _sum_sq_charge << std::endl;
      std::cout << "_qqr2e: " << _qqr2e << std::endl;
      std::cout << "_q2: " << _q2 << std::endl;

      //compute g_ewald
      _g_ewald = _accuracy*SQRT(num_atoms * _cut_off * volume) / (2.0 * _q2);
      if (_g_ewald >= 1.0) _g_ewald = (1.35 - 0.15*LOG(_accuracy))/_cut_off;
      else _g_ewald = SQRT(-LOG(_g_ewald)) / _cut_off;
      _alpha = _g_ewald*_g_ewald;
      std::cout << "relative_accuracy_user : " <<relative_accuracy_user<<
        ", _g_ewald; "<<_g_ewald  << ",  alpha= " <<  _alpha <<std::endl;


    }

    //alpha
    if (config->PathExists({"hyper_parameters", "coulomb" ,"alpha"})) {
      _alpha = config->Get<rbmd::Real>("alpha", "hyper_parameters", "coulomb");

      auto accuracy_test = ERFC(_cut_off * SQRT(_alpha));
      // std::cout << "accuracy_test= " <<  accuracy_test <<std::endl;
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
    //RBSOG
    else if ("RBSOG" == _coulomb_type)
    {
      // Read RBSOG specific parameters
      _RBE_P = config->Get<rbmd::Id>("coulomb_sample_num", "hyper_parameters", "coulomb"); // This is 'P'
      _rbsog_b = config->Get<rbmd::Real>("rbsog_b", "hyper_parameters", "coulomb");
      _rbsog_sigma = config->Get<rbmd::Real>("rbsog_sigma", "hyper_parameters", "coulomb");
      _rbsog_Mmax = config->Get<rbmd::Id>("rbsog_Mmax", "hyper_parameters", "coulomb");
      _rbsog_Kcut = config->Get<rbmd::Id>("rbsog_Kcut", "hyper_parameters", "coulomb");
      _rbsog_Kcut = 1;

      _h_rbsog_K_Sample_int.resize(_RBE_P);
      _idx_npt.resize(_RBE_P,0);
      _h_rbsog_idx_npt.resize(_RBE_P,0);

      // Initialize RBSOG coefficients and S-values
      RBSOGInit();
      RBSOGSetup(); // Initial calculation of S, S_npt
    }
  }
}

void KSpaceCalculator::Execute()
{
    if ("RBE" == _coulomb_type)
    {
        ComputeRBE();
    }
    else if ("RBSOG" == _coulomb_type)
    {
        ComputeRBSOG();
    }
    else //
    {
        ComputeEwald();
    }

    EvaluatePotentialEnergy();
}

void KSpaceCalculator::EvaluatePotentialEnergy()
{
  auto num_atoms = *(_structure_info_data->_num_atoms);

  auto unit = DataManager::getInstance().getConfigData()->Get
<std::string>("unit", "init_configuration", "read_data");
  if ("LJ" == unit) {
    _e_kspace = _e_kspace/num_atoms;
  }
}

void KSpaceCalculator::ComputeEwald()
{
    auto start = std::chrono::high_resolution_clock::now();
    auto num_atoms = *(_structure_info_data->_num_atoms);

    //compute charge structure factor for Ewald
    thrust::host_vector<rbmd::Real> h_Re_array(_num_k);
    thrust::host_vector<rbmd::Real> h_Im_array(_num_k);
    ComputeChargeStructureFactorEwald(*_box, num_atoms, _kmax_array, _alpha,
      _qqr2e, h_Re_array, h_Im_array);

    // H2D
    thrust::device_vector<rbmd::Real> d_real_array = h_Re_array;
    thrust::device_vector<rbmd::Real> d_imag_array = h_Im_array;

    //EwaldForce on GPU//
    op::ComputeEwaldForceOp<device::DEVICE_GPU>()(
        *_box, num_atoms, _kmax_array, _alpha, _qqr2e,
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
    TimingStatistics::Instance().record("Long-Range", duration.count());

  // thrust::host_vector<rbmd::Real> h_kspace_x =_device_data->_d_force_kspace_x;
  // thrust::host_vector<rbmd::Real> h_kspace_y = _device_data->_d_force_kspace_y;
  // thrust::host_vector<rbmd::Real> h_kspace_z = _device_data->_d_force_kspace_z;
  //
  //
  //
  // std::ofstream kspace_file("kspace_ewald.txt");
  // auto atom_id_to_idx =
  //   LinkedCellLocator::GetInstance().GetLinkedCell()->_atom_id_to_idx;
  // if (kspace_file.is_open()) {
  //   for (rbmd::Id i = 0; i < h_kspace_x.size(); ++i) {
  //     auto index = atom_id_to_idx[i];
  //     kspace_file << i  << " " <<h_kspace_x[index] <<" " <<h_kspace_y[index] <<" " <<
  //     h_kspace_z[index]<< "\n";
  //   }
  //   kspace_file.close();
  // }

    //
    ReduceVirial(num_atoms, _device_data->_d_flat_virial_kspace,
      _device_data->_d_virial_kspace);

  // thrust::host_vector<rbmd::Real> h_kspace_ewald_virial =_device_data->_d_virial_kspace;
  //
  // std::ofstream kspace_ewald_file("kspace_ewald_virial.txt");
  // if (kspace_ewald_file.is_open()) {
  //   for (rbmd::Id i = 0; i < h_kspace_ewald_virial.size(); ++i) {
  //     kspace_ewald_file << i  << " " <<h_kspace_ewald_virial[i]  << "\n";
  //   }
  //   kspace_ewald_file.close();
  // }
}

void KSpaceCalculator::ComputeChargeStructureFactorEwald(
    Box box, rbmd::Id num_atoms, Int3 kmax_array,
    rbmd::Real alpha, rbmd::Real qqr2e,
    thrust::host_vector<rbmd::Real>& value_Re_array,
    thrust::host_vector<rbmd::Real>& value_Im_array)
{
    // Performance tip:  Launching a CUDA kernel for each k vector in the triple loop and performing reduce is inefficient.
    // A more optimized approach is to write a single CUDA kernel to process all k vectors in parallel.

    thrust::device_vector<rbmd::Real> density_real_atom(num_atoms);
    thrust::device_vector<rbmd::Real> density_imag_atom(num_atoms);

    rbmd::Real total_energy_kspace = 0;
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
                    rbmd::Real Range_K2 = K.x * K.x + K.y * K.y + K.z * K.z;
                    rbmd::Real alpha_inv = 1.0 / alpha;

                    op::ComputeChargeStructureFactorOp<device::DEVICE_GPU>()(
                        num_atoms, K,
                        thrust::raw_pointer_cast(_device_data->_d_charge.data()),
                        thrust::raw_pointer_cast(_device_data->_d_px.data()),
                        thrust::raw_pointer_cast(_device_data->_d_py.data()),
                        thrust::raw_pointer_cast(_device_data->_d_pz.data()),
                        thrust::raw_pointer_cast(density_real_atom.data()),
                        thrust::raw_pointer_cast(density_imag_atom.data()));

                    rbmd::Real value_Re = thrust::reduce(density_real_atom.begin(), density_real_atom.end());
                    rbmd::Real value_Im = thrust::reduce(density_imag_atom.begin(), density_imag_atom.end());
                    rbmd::Real Range_density2 = value_Re * value_Re + value_Im * value_Im;

                    total_energy_kspace += EXP(-0.25 * Range_K2 * alpha_inv) * Range_density2 / Range_K2;

                    value_Re_array[index] = value_Re;
                    value_Im_array[index] = value_Im;
                    index++;
                }
            }
        }
    }

   //energy
   //charge self energy//
    ComputeSelfEnergy(alpha, qqr2e, _e_self_energy);

   //compute Kspace energy//
    rbmd::Real volume = box._length[0] * box._length[1] * box._length[2];
    total_energy_kspace = qqr2e * (2 * M_PI / volume) * total_energy_kspace;
    _e_kspace = total_energy_kspace  + _e_self_energy;
}


void KSpaceCalculator::ComputeRBE()
{
  auto start = std::chrono::high_resolution_clock::now();
  //
  auto num_atoms = *(_structure_info_data->_num_atoms);
  _rhok_real_redue.resize(_RBE_P);
  _rhok_image_redue.resize(_RBE_P);

  ComputeChargeStructureFactorRBE_opt(*_box, num_atoms, _kmax_array,
      _alpha,_RBE_P,_qqr2e,_rhok_real_redue,_rhok_image_redue);

  //RBE Force
  op::ComputeRBEForceOp<device::DEVICE_GPU>()(
        *_box,num_atoms, _RBE_P,_alpha,_qqr2e,
        thrust::raw_pointer_cast(_rhok_real_redue.data()),
        thrust::raw_pointer_cast(_rhok_image_redue.data()),
        thrust::raw_pointer_cast(_device_data->_d_charge.data()),
        thrust::raw_pointer_cast(_P_Sample_x.data()),
        thrust::raw_pointer_cast(_P_Sample_y.data()),
        thrust::raw_pointer_cast(_P_Sample_z.data()),
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


  auto start1 = std::chrono::high_resolution_clock::now();
  ComputeRBEVirial();
  auto end1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration1 = end1 - start1;
  // std::cout << "time of old ComputeRBEVirial on GPU: "<<duration1.count() << std::endl;

  //   ReduceVirial(num_atoms,_device_data->_d_flat_virial_kspace,
  // _device_data->_d_virial_kspace);

  // thrust::host_vector<rbmd::Real> h_kspace_x = _device_data->_d_force_kspace_x;
  // thrust::host_vector<rbmd::Real> h_kspace_y = _device_data->_d_force_kspace_y;
  // thrust::host_vector<rbmd::Real> h_kspace_z= _device_data->_d_force_kspace_z;
  // std::ofstream kspace_file("kspace_rbe.txt");
  // auto atom_id_to_idx =
  //   LinkedCellLocator::GetInstance().GetLinkedCell()->_atom_id_to_idx;
  // if (kspace_file.is_open()) {
  //   for (rbmd::Id i = 0; i < h_kspace_x.size(); ++i) {
  //     auto index = atom_id_to_idx[i];
  //     kspace_file << i  << " " <<h_kspace_x[index] <<" " <<h_kspace_y[index] <<" " <<
  //     h_kspace_z[index]<< "\n";
  //   }
  //   kspace_file.close();
  // }


}

void KSpaceCalculator::ComputeChargeStructureFactorRBE(
   Box box,
   rbmd::Id num_atoms,
   Int3 kmax_array,
   rbmd::Real alpha,
   rbmd::Id RBE_P,
   rbmd::Real qqr2e,
   thrust::device_vector<rbmd::Real>& rhok_real_redue,
   thrust::device_vector<rbmd::Real>& rhok_image_redue)
{
  auto start = std::chrono::high_resolution_clock::now();
  //get P_Sample at each step
  RBEInit(*_box,_alpha,_RBE_P);

  thrust::device_vector<rbmd::Real>  rhok_real_atom;
  thrust::device_vector<rbmd::Real>  rhok_image_atom;
  rhok_real_atom.resize(num_atoms* RBE_P);
  rhok_image_atom.resize(num_atoms* RBE_P);
  auto p_number= RBE_P;

  //Charge Structure Factor for RBE
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

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration = end - start;
  std::cout <<  "传统结构因子耗时："<< duration.count() << std::endl;

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

void KSpaceCalculator::ComputeChargeStructureFactorRBE_opt(
   Box box,
   rbmd::Id num_atoms,
   Int3 kmax_array,
   rbmd::Real alpha,
   rbmd::Id RBE_P,
   rbmd::Real qqr2e,
   thrust::device_vector<rbmd::Real>& rhok_real_redue,
   thrust::device_vector<rbmd::Real>& rhok_image_redue)
{
  auto start = std::chrono::high_resolution_clock::now();
  //get P_Sample at each step
  RBEInit(*_box,_alpha,_RBE_P);
  auto p_number= RBE_P;

  //Charge Structure Factor for RBE
  op::ComputePnumberChargeStructureFactorOp<device::DEVICE_GPU>()(
      box, num_atoms, p_number,
      thrust::raw_pointer_cast(_device_data->_d_charge.data()),
      thrust::raw_pointer_cast(_P_Sample_x.data()),
      raw_pointer_cast(_P_Sample_y.data()),
      raw_pointer_cast(_P_Sample_z.data()),
      thrust::raw_pointer_cast(_device_data->_d_px.data()),
      thrust::raw_pointer_cast(_device_data->_d_py.data()),
      thrust::raw_pointer_cast(_device_data->_d_pz.data()),
      thrust::raw_pointer_cast(rhok_real_redue.data()),
      thrust::raw_pointer_cast(rhok_image_redue.data()));

  // auto end = std::chrono::high_resolution_clock::now();
  // std::chrono::duration<rbmd::Real> duration = end - start;
  // std::cout <<  "优化结构因子耗时："<< duration.count() << std::endl;


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

void KSpaceCalculator::ComputeRBEVirial()
{
  //
  auto num_atoms = *(_structure_info_data->_num_atoms);
  rbmd::Real h_energy;
  thrust::device_vector<rbmd::Real> d_energy(1);

  //
  op::ComputeRBEForceVirialOp<device::DEVICE_GPU>()(
        *_box ,_RBE_P,_alpha,_qqr2e,_sum_sq_charge,_sum_charge,
        thrust::raw_pointer_cast(_rhok_real_redue.data()),
        thrust::raw_pointer_cast(_rhok_image_redue.data()),
        thrust::raw_pointer_cast(_P_Sample_x.data()),
        thrust::raw_pointer_cast(_P_Sample_y.data()),
        thrust::raw_pointer_cast(_P_Sample_z.data()),
        thrust::raw_pointer_cast(_device_data->_d_virial_kspace.data()),
        thrust::raw_pointer_cast( d_energy.data()));

  // D2H
  thrust::copy(d_energy.begin(), d_energy.end(), &h_energy);
  _e_kspace = h_energy;

  // thrust::host_vector<rbmd::Real> h_kspace_rbe_virial =_device_data->_d_virial_kspace;
  //
  // std::ofstream kspace_rbe_file("kspace_rbe_virial.txt");
  // if (kspace_rbe_file.is_open()) {
  //   for (rbmd::Id i = 0; i < h_kspace_rbe_virial.size(); ++i) {
  //     kspace_rbe_file << i  << " " <<h_kspace_rbe_virial[i]  << "\n";
  //   }
  //   kspace_rbe_file.close();
  // }

}

void KSpaceCalculator::RBEInit(Box box,rbmd::Real alpha,rbmd::Id RBE_P)
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

void KSpaceCalculator::GetPsampleKey()
{
  //psample index key
  auto num_atoms = *(_structure_info_data->_num_atoms);
  _psample_key.resize(num_atoms * _RBE_P);
  op::GenerateIndexArrayOp<device::DEVICE_GPU>()(num_atoms,_RBE_P,
    thrust::raw_pointer_cast(_psample_key.data()));
}

void KSpaceCalculator::ComputeQsqSum()
{
  auto num_atoms = *(_structure_info_data->_num_atoms);
  thrust::device_vector<rbmd::Real> sq_charge;
  sq_charge.resize(num_atoms);
  op::SqchargeOp<device::DEVICE_GPU>()(num_atoms,
    thrust::raw_pointer_cast(_device_data->_d_charge.data()),
    thrust::raw_pointer_cast(sq_charge.data()));

  _sum_sq_charge = thrust::reduce(sq_charge.begin(),
   sq_charge.end(), 0.0f, thrust::plus<rbmd::Real>());
  _q2 = _qqr2e* _sum_sq_charge;

  //
  _sum_charge = thrust::reduce(_device_data->_d_charge.begin(),
_device_data->_d_charge.end(), 0.0f, thrust::plus<rbmd::Real>());

  if (ABS(_sum_charge) > SMALL) {
    Logger::Instance().warn("\033[31mThe total charge of the model is not zero, "
              "the net charge is {}.\033[0m", _sum_charge);
  }
}

void KSpaceCalculator::ComputeKspaceEnergy(
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
  ave_ekspace = total_energy_ewald ;
}

void KSpaceCalculator::ComputeSelfEnergy(
  rbmd::Real  alpha,
  rbmd::Real  qqr2e,
  rbmd::Real& ave_self_energy)
{
  //compute self_energy
  auto num_atoms = *(_structure_info_data->_num_atoms);
  rbmd::Real total_self_energy = qqr2e * (- SQRT(alpha / M_PI) *_sum_sq_charge);

  ave_self_energy =  total_self_energy ;
}

rbmd::Real KSpaceCalculator::ComputeRMS(rbmd::Id kmax,rbmd::Real box_length,rbmd::Real q2)
{
  auto num_atoms = *(_structure_info_data->_num_atoms);
  rbmd::Real rms = 2.0*q2*_g_ewald/box_length *
    SQRT(1.0/(M_PI*kmax*num_atoms)) *
    EXP(-M_PI*M_PI*kmax*kmax/(_g_ewald*_g_ewald*box_length*box_length));

  return rms;
}

void KSpaceCalculator::SetKspacePara()
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

  rms = ComputeRMS(kmax_x,box->_length[0],_q2);
  while (rms > _accuracy) {
    kmax_x++;
    rms = ComputeRMS(kmax_x,box->_length[0],_q2);
  }

  rms = ComputeRMS(kmax_y,box->_length[1],_q2);
  while (rms > _accuracy) {
    kmax_y++;
    rms = ComputeRMS(kmax_y,box->_length[1],_q2);
  }

  rms = ComputeRMS(kmax_z,box->_length[2],_q2);
  while (rms > _accuracy) {
    kmax_z++;
    rms = ComputeRMS(kmax_z,box->_length[2],_q2);
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

void KSpaceCalculator::RBSOGInit()
{
  rbmd::Real r0 = _cut_off / _rbsog_sigma; //归一化截断半径
  _rbsog_w0 = Compute_W0(r0, _rbsog_b); // 计算最窄高斯权重（式2.11）   b: 几何级数的基，控制高斯带宽的间隔

  const auto& config = DataManager::getInstance().getConfigData();
  if (config->PathExists({"hyper_parameters", "coulomb" ,"rbsog_omega"})) {
    _rbsog_omega = config->Get<rbmd::Real>("rbsog_omega", "hyper_parameters", "coulomb");
    _rbsog_w0= _rbsog_omega;
  }
  std::cout <<  "_rbsog_w0: " << _rbsog_w0  <<std::endl;

  _h_rbsog_sl.resize(_rbsog_Mmax);
  _h_rbsog_coef.resize(_rbsog_Mmax);
  _h_rbsog_coef_npt.resize(_rbsog_Mmax);

  _h_rbsog_sl[0] = SQRT(2.0) * _rbsog_sigma; // 高斯宽度s_l（式2.9） ,乘以√2以匹配文献中的归一化形式。当l=0时。。。sl[0] = √2* sigma
  rbmd::Real sigma2 = _rbsog_sigma * _rbsog_sigma;
  rbmd::Real sigma4 = sigma2 * sigma2;
  rbmd::Real b2 = _rbsog_b * _rbsog_b;
  rbmd::Real b4 = b2 * b2;

  _h_rbsog_coef[0] = 4 * M_PI * LOG(_rbsog_b) * _rbsog_w0 * sigma2; // 权重系数初始化, 第一个高斯项 3.5 在远场势计算中的权重系数，
  _h_rbsog_coef_npt[0] = 8 * M_PI * LOG(_rbsog_b) * _rbsog_w0 * sigma4;

  // 循环计算各高斯项参数（式2.8-2.9）
  for (int i = 1; i < _rbsog_Mmax; ++i) {
      _h_rbsog_sl[i] = _h_rbsog_sl[i-1] * _rbsog_b;   //（递推几何级数）
      // This logic from rbsog_intel.cpp seems slightly off (i==1 block looks redundant)
      // but porting it faithfully:
      if( i == 1){
          _h_rbsog_coef[i] = 4 * M_PI * LOG(_rbsog_b) * sigma2 * b2;
          _h_rbsog_coef_npt[i] = 8 * M_PI * LOG(_rbsog_b) * sigma4 * b4;
      }
      _h_rbsog_coef[i] = _h_rbsog_coef[i-1] * b2;
      _h_rbsog_coef_npt[i] = _h_rbsog_coef_npt[i-1] * b4;
  }

  // 1.
  // std::ofstream sl_file("_h_rbsog_sl.txt");
  // if (sl_file.is_open()) {
  //   for (rbmd::Id i = 0; i < _h_rbsog_sl.size(); ++i) {
  //     sl_file << i  << " " << _h_rbsog_sl[i] << "\n";
  //   }
  //   sl_file.close();
  // }
  //
  // std::ofstream coef_file("_h_rbsog_coef.txt");
  // if (coef_file.is_open()) {
  //   for (rbmd::Id i = 0; i < _h_rbsog_coef.size(); ++i) {
  //     coef_file << i  << " " << _h_rbsog_coef[i] << "\n";
  //   }
  //   coef_file.close();
  // }

  // Copy coefficients to device for kernel use
  _d_rbsog_sl = _h_rbsog_sl;
  _d_rbsog_coef = _h_rbsog_coef;
  _d_rbsog_coef_npt = _h_rbsog_coef_npt;
}

void KSpaceCalculator::RBSOGSetup()
{
  rbmd::Real xprd = _box->_length[0];
  rbmd::Real yprd = _box->_length[1];
  rbmd::Real zprd = _box->_length[2];
  rbmd::Real prd = MAX(xprd, MAX(yprd, zprd));

  rbmd::Real sum = 0.0;
  rbmd::Real sum_npt = 0.0;
  rbmd::Real sum_1 = 0.0;
  rbmd::Real sum_npt_1 = 0.0;

  rbmd::Real sigma2 = _rbsog_sigma * _rbsog_sigma;
  rbmd::Real sigma4 = sigma2 * sigma2;
  rbmd::Real log_b = LOG(_rbsog_b);

  // Direct sum part (low-k)
  rbmd::Real bool_bound = 1.0;
  for (int i = 0; i < _rbsog_Kcut +1 ; i++)
  {
        rbmd::Real Kx = (2 * M_PI / xprd) * (i + 0.00);
        for (int j = 0; j < _rbsog_Kcut + 1; j++)
        {
            rbmd::Real Ky = (2 * M_PI / yprd) * (j + 0.00);
            for (int k = 0; k < _rbsog_Kcut + 1; k++)
            {
                rbmd::Real Kz = (2 * M_PI / zprd) * (k + 0.00);
                if ((i * i + j * j + k * k) <= _rbsog_Kcut)
                {
                    rbmd::Real K2 = Kx * Kx + Ky * Ky + Kz * Kz;
                    rbmd::Real K4 = K2 * K2;
                    sum_1 += 2*Gaussian(i, j, k, *_box, _rbsog_sigma, _rbsog_b, _rbsog_w0, _rbsog_Mmax, _h_rbsog_coef) * K2; //低波矢离散求和
                    sum_npt_1 += 2*Gaussian_modify(i, j, k,  *_box, _rbsog_sigma, _rbsog_b, _rbsog_w0, _rbsog_Mmax, _h_rbsog_coef_npt) * K4;
                }
            }
        }
  }


  // Analytical integral part (high-k)
  for (int i = 0; i < _rbsog_Mmax; ++i)
  {
      rbmd::Real Hx = 0, Hy = 0, Hz = 0; // 累积高斯项
      rbmd::Real Yx = 0, Yy = 0, Yz = 0; // 累积一阶导数项
      rbmd::Real Yxx = 0, Yyy = 0, Yzz = 0; // 累积二阶导数项

      int j = 1;  //波矢索引（从1开始）
      while (bool_bound > 1e-12) {
          rbmd::Real kx = M_PI * j / xprd;
          rbmd::Real ky = M_PI * j / yprd;
          rbmd::Real kz = M_PI * j / zprd;
          rbmd::Real k_min = M_PI * j / prd;
          rbmd::Real k_min2 = k_min * k_min;
          rbmd::Real k_min4 = k_min2 * k_min2;

          rbmd::Real sl_i = _h_rbsog_sl[i];
          rbmd::Real sl_i_sq = sl_i * sl_i;

          rbmd::Real kx2_s = kx * sl_i * kx * sl_i; // This is (kx*sl[i])^2
          rbmd::Real ky2_s = ky * sl_i * ky * sl_i;
          rbmd::Real kz2_s = kz * sl_i * kz * sl_i;

          rbmd::Real kx2 = 4 * kx * kx; // This is (2*kx)^2
          rbmd::Real ky2 = 4 * ky * ky;
          rbmd::Real kz2 = 4 * kz * kz;
          rbmd::Real kx4 = 16 * kx * kx * kx * kx; // This is (2*kx)^4
          rbmd::Real ky4 = 16 * ky * ky * ky * ky;
          rbmd::Real kz4 = 16 * kz * kz * kz * kz;

          Hx += 2 * EXP(-kx2_s);
          Yx += 2 * EXP(-kx2_s) * kx2;
          Yxx += 2 * EXP(-kx2_s) * kx4;
          Hy += 2 * EXP(-ky2_s);
          Yy += 2 * EXP(-ky2_s) * ky2;
          Yyy += 2 * EXP(-ky2_s) * ky4;
          Hz += 2 * EXP(-kz2_s);
          Yz += 2 * EXP(-kz2_s) * kz2;
          Yzz += 2 * EXP(-kz2_s) * kz4;
          ++j;

          // Original bool_bound logic
          bool_bound = EXP(4 * i * log_b - k_min2 * sl_i_sq) * k_min4 * sigma4;
      }
      Hx += 1; // Add j=0 term
      Hy += 1;
      Hz += 1;

      sum += _h_rbsog_coef[i] * (Yx * Hy * Hz + Hx * Yy * Hz + Hx * Hy * Yz); // 累加总和高波矢贡献,高波矢解析积分
      sum_npt += _h_rbsog_coef_npt[i] * (Yxx * Hy * Hz + Hx * Yyy * Hz + Hx * Hy * Yzz + 2 * Yx * Yy * Hz + 2 * Yx * Hy * Yz + 2 * Hx * Yy * Yz);
  }

  _rbsog_S = sum - sum_1; //消除重复计算部分，确保全波矢空间的一致性。
  _rbsog_S_npt = sum_npt - sum_npt_1;

  std::cout<< "_rbsog_S:" << _rbsog_S << ",_rbsog_S_npt:" << _rbsog_S_npt<<std::endl;
  // Logger::Instance().info("RBSOGSetup: S = " + std::to_string(_rbsog_S) + ", S_npt = " + std::to_string(_rbsog_S_npt));
}

void KSpaceCalculator::RBSOGSampleKSpace()
{
  _h_rbsog_K_Sample_x.resize(_RBE_P);
  _h_rbsog_K_Sample_y.resize(_RBE_P);
  _h_rbsog_K_Sample_z.resize(_RBE_P);

  rbmd::Real xprd = _box->_length[0];
  rbmd::Real yprd = _box->_length[1];
  rbmd::Real zprd = _box->_length[2];

  rbmd::Real pxyz[3] = { (rbmd::Real)(2 * M_PI / xprd),
                        (rbmd::Real)(2 * M_PI / yprd),
                        (rbmd::Real)(2 * M_PI / zprd) };

  // Host vectors for sampling
  thrust::host_vector<rbmd::Id> mx,my,mz;
  mx.resize(_RBE_P * 2);
  my.resize(_RBE_P * 2);
  mz.resize(_RBE_P * 2);


  rbmd::Real factor_xyz[3] = { rbmd::Real(xprd / (2 * M_PI * _rbsog_sigma)) ,
                               rbmd::Real(yprd / (2 * M_PI * _rbsog_sigma)) ,
                               rbmd::Real(zprd / (2 * M_PI * _rbsog_sigma))};
  rbmd::Real MHD_factor[3] = { SQRT(2 * _rbsog_sigma * _rbsog_sigma * M_PI * M_PI) / xprd,
                               SQRT(2 * _rbsog_sigma * _rbsog_sigma * M_PI * M_PI) / yprd,
                               SQRT(2 * _rbsog_sigma * _rbsog_sigma * M_PI * M_PI) / zprd };

  rbmd::Real x, mold_x, mnew_x, xx, mold_y, mnew_y, xxx, mold_z, mnew_z;
  rbmd::Real Kx_new, Ky_new, Kz_new, K2_new, Kx_old, Ky_old, Kz_old, K2_old;
  rbmd::Real pup, qup, pdown, qdown, acce, yyy;

  // 1. Generate MCMC chain for K_Sample (Force)
  do {
      mx[0] = ROUND(randn_box_muller(0, factor_xyz[0]));
      my[0] = ROUND(randn_box_muller(0, factor_xyz[1]));
      mz[0] = ROUND(randn_box_muller(0, factor_xyz[2]));
  } while (mx[0] == 0 && my[0] == 0 && mz[0] == 0);

  for (int i = 0; i < 2 * _RBE_P - 1; i++)
  {
      x = randn_box_muller(0, factor_xyz[0]);
      mold_x = mx[i];
      mnew_x = ROUND(x);

      xx = randn_box_muller(0, factor_xyz[1]);
      mold_y = my[i];
      mnew_y = ROUND(xx);

      xxx = randn_box_muller(0, factor_xyz[2]);
      mold_z = mz[i];
      mnew_z = ROUND(xxx);

      Kx_new = pxyz[0] * (mnew_x + 0.00);
      Ky_new = pxyz[1] * (mnew_y + 0.00);
      Kz_new = pxyz[2] * (mnew_z + 0.00);
      K2_new = Kx_new * Kx_new + Ky_new * Ky_new + Kz_new * Kz_new;
      Kx_old = pxyz[0] * (mold_x + 0.00);
      Ky_old = pxyz[1] * (mold_y + 0.00);
      Kz_old = pxyz[2] * (mold_z + 0.00);
      K2_old = Kx_old * Kx_old + Ky_old * Ky_old + Kz_old * Kz_old;

      pup = Gaussian((int)mnew_x, (int)mnew_y, (int)mnew_z,  *_box, _rbsog_sigma,
        _rbsog_b, _rbsog_w0, _rbsog_Mmax, _h_rbsog_coef) * K2_new;
      qup = MH_D_Modify((int)mold_x, MHD_factor[0]) * MH_D_Modify((int)mold_y,
        MHD_factor[1]) * MH_D_Modify((int)mold_z, MHD_factor[2]);

      pdown = Gaussian((int)mold_x, (int)mold_y, (int)mold_z,  *_box,
        _rbsog_sigma, _rbsog_b, _rbsog_w0, _rbsog_Mmax, _h_rbsog_coef) * K2_old;
      qdown = MH_D_Modify((int)mnew_x, MHD_factor[0]) * MH_D_Modify((int)mnew_y,
        MHD_factor[1]) * MH_D_Modify((int)mnew_z, MHD_factor[2]);

      acce = pup * qup / (pdown * qdown) > 1.0 ? 1.0 : pup * qup / (pdown * qdown);

      yyy = RandomValue< rbmd::Real>(0.0, 1.0);

      if (yyy < acce) {
          mx[i + 1] = mnew_x;
          my[i + 1] = mnew_y;
          mz[i + 1] = mnew_z;
      }
      else {
          mx[i + 1] = mold_x;
          my[i + 1] = mold_y;
          mz[i + 1] = mold_z;
      }

      int m_sq = mx[i+1]*mx[i+1] + my[i+1]*my[i+1] + mz[i+1]*mz[i+1];
      if ((mx[i + 1]  == 0 && my[i + 1] == 0 && mz[i + 1] == 0) || (m_sq <= _rbsog_Kcut))
      {
          i = i - 1; // Reject sample and retry
      }
  }

  // 2. Generate MCMC chain for idx_npt (Virial)
  // This chain depends on the K_Sample chain

  for (int i = 0; i < _RBE_P; i++)
  {
    _h_rbsog_K_Sample_x[i]  = mx[2 * i + 1] + 0.00;
    _h_rbsog_K_Sample_y[i]  = my[2 * i + 1] + 0.00;
    _h_rbsog_K_Sample_z[i]  = mz[2 * i + 1] + 0.00;
  }

  for (int i = 1; i < _RBE_P; i++)
  {
      rbmd::Real mmx_old, mmy_old, mmz_old, mmx_new, mmy_new, mmz_new;
      rbmd::Real K4_new, K2_old, K4_old, K2_new;
      rbmd::Real pup, qup, pdown, qdown, acce, yyy;

      rbmd::Id idx;
      idx = _idx_npt[i - 1]; //

      // Get K_old from K_Sample chain at previous accepted index
      mmx_old =  _h_rbsog_K_Sample_x[idx] ;
      mmy_old =  _h_rbsog_K_Sample_y[idx] ;
      mmz_old =  _h_rbsog_K_Sample_z[idx] ;

      // Get K_new from K_Sample chain at current index
      mmx_new =  _h_rbsog_K_Sample_x[i] ;
      mmy_new =  _h_rbsog_K_Sample_y[i] ;
      mmz_new =  _h_rbsog_K_Sample_z[i] ;

      Kx_new = pxyz[0] * mmx_new;
      Ky_new = pxyz[1] * mmy_new;
      Kz_new = pxyz[2] * mmz_new;
      K2_new = Kx_new * Kx_new + Ky_new * Ky_new + Kz_new * Kz_new;
      K4_new = K2_new * K2_new;
      Kx_old = pxyz[0] * mmx_old;
      Ky_old = pxyz[1] * mmy_old;
      Kz_old = pxyz[2] * mmz_old;
      K2_old = Kx_old * Kx_old + Ky_old * Ky_old + Kz_old * Kz_old;
      K4_old = K2_old * K2_old;

      pup = Gaussian_modify((int)mmx_new, (int)mmy_new, (int)mmz_new,  *_box,
        _rbsog_sigma, _rbsog_b, _rbsog_w0, _rbsog_Mmax, _h_rbsog_coef_npt) * K4_new;
      qup = Gaussian((int)mmx_old, (int)mmy_old, (int)mmz_old,  *_box,
        _rbsog_sigma, _rbsog_b, _rbsog_w0, _rbsog_Mmax, _h_rbsog_coef) * K2_old;

      pdown = Gaussian_modify((int)mmx_old, (int)mmy_old, (int)mmz_old,
         *_box, _rbsog_sigma, _rbsog_b, _rbsog_w0, _rbsog_Mmax, _h_rbsog_coef_npt) * K4_old;
      qdown = Gaussian((int)mmx_new, (int)mmy_new, (int)mmz_new,
         *_box, _rbsog_sigma, _rbsog_b, _rbsog_w0, _rbsog_Mmax, _h_rbsog_coef) * K2_new;

      acce = pup * qup / (pdown * qdown) > 1.0 ? 1.0 : pup * qup / (pdown * qdown);

      yyy = RandomValue< rbmd::Real>(0.0, 1.0);

      if (yyy < acce) {
          _idx_npt[i] = i; // Accept current sample index
      }
      else {
          _idx_npt[i] = idx; // Reject, reuse previous sample index
      }
  }

  // 3. Prepare device vectors
  // We need to store the *real* K-vectors, not the integer (m) vectors
  thrust::host_vector<rbmd::Real> h_K_x(_RBE_P), h_K_y(_RBE_P), h_K_z(_RBE_P);
  thrust::host_vector<rbmd::Real> h_K_npt_x(_RBE_P), h_K_npt_y(_RBE_P), h_K_npt_z(_RBE_P);

  for (int i = 0; i < _RBE_P; i++)
  {
      // K_Sample (force) vectors (using the "burned-in" samples)
      h_K_x[i] =_h_rbsog_K_Sample_x[i] * pxyz[0];
      h_K_y[i] =_h_rbsog_K_Sample_y[i] * pxyz[1];
      h_K_z[i] =_h_rbsog_K_Sample_z[i] * pxyz[2];

      // K_npt (virial) vectors (gathered from the force list)
      rbmd::Id id = _idx_npt[i];
      h_K_npt_x[i] = _h_rbsog_K_Sample_x[id] * pxyz[0];
      h_K_npt_y[i] = _h_rbsog_K_Sample_y[id] * pxyz[1];
      h_K_npt_z[i] = _h_rbsog_K_Sample_z[id] * pxyz[2];
  }

  // std::ofstream K_Sample_file("h_rbsog_K_Sample_int.txt");
  // if (K_Sample_file.is_open()) {
  //   for (rbmd::Id i = 0; i < _RBE_P; ++i) {
  //     K_Sample_file << i  << " " << _h_rbsog_K_Sample_int[i].x<< " "
  //     <<  _h_rbsog_K_Sample_int[i].y<< " "
  //     <<  _h_rbsog_K_Sample_int[i].z<<"\n";
  //   }
  //   K_Sample_file.close();
  // }


  // 4. Copy to device
  _d_rbsog_K_Sample_x = h_K_x;
  _d_rbsog_K_Sample_y = h_K_y;
  _d_rbsog_K_Sample_z = h_K_z;
  _d_rbsog_K_npt_x = h_K_npt_x;
  _d_rbsog_K_npt_y = h_K_npt_y;
  _d_rbsog_K_npt_z = h_K_npt_z;
  _d_rbsog_idx_npt_all = _idx_npt; // This is the gather map
}

void KSpaceCalculator::ComputeRBSOG() {
  auto start = std::chrono::high_resolution_clock::now();
  auto num_atoms = *(_structure_info_data->_num_atoms);
  rbmd::Real xprd = _box->_length[0];
  rbmd::Real yprd = _box->_length[1];
  rbmd::Real zprd = _box->_length[2];

  // // if(currt_step % interl  ==0 )
  xprd0 = xprd;
  yprd0 = yprd;
  zprd0 = zprd;
  _rbsog_S0 = _rbsog_S;
  _rbsog_S_npt0 = _rbsog_S_npt;

  // 2. Perform MCMC sampling (host-side) and upload K-vectors to device
  RBSOGSampleKSpace();  //

  // 3. Calculate correction factors `fac` and `fac_npt` on GPU
  // auto start1 = std::chrono::high_resolution_clock::now();
  _d_fac.resize(_RBE_P);
  _d_fac_npt.resize(_RBE_P);
  thrust::fill(_d_fac.begin(), _d_fac.end(), 1);
  thrust::fill(_d_fac_npt.begin(), _d_fac_npt.end(), 1);

  //
  rbmd::Real L_ratio = xprd / xprd0;
  rbmd::Real S_ratio = _rbsog_S0 / _rbsog_S;
  rbmd::Real S_npt_ratio = _rbsog_S_npt0 / _rbsog_S_npt;
  // ComputeRBSOGFactor();
  op::ComputeRBSOGFactorsOp<device::DEVICE_GPU>()(
      _RBE_P, thrust::raw_pointer_cast(_d_rbsog_K_Sample_x.data()),
      thrust::raw_pointer_cast(_d_rbsog_K_Sample_y.data()),
      thrust::raw_pointer_cast(_d_rbsog_K_Sample_z.data()),
      thrust::raw_pointer_cast(_d_rbsog_K_npt_x.data()),
      thrust::raw_pointer_cast(_d_rbsog_K_npt_y.data()),
      thrust::raw_pointer_cast(_d_rbsog_K_npt_z.data()), *_box, _rbsog_sigma,
      _rbsog_b, _rbsog_Mmax, thrust::raw_pointer_cast(_d_rbsog_coef.data()),
      thrust::raw_pointer_cast(_d_rbsog_coef_npt.data()), L_ratio, S_ratio,
      S_npt_ratio, thrust::raw_pointer_cast(_d_fac.data()),
      thrust::raw_pointer_cast(_d_fac_npt.data()));
  // auto end1 = std::chrono::high_resolution_clock::now();
  // std::chrono::duration<rbmd::Real> duration1 = end1 - start1;
  // std::cout << "RBSOGFactor time : " << duration1.count() << std::endl;

  // thrust::host_vector<rbmd::Real> h_fac = _d_fac;
  // std::ofstream fac_file("h_fac.txt");
  // if (fac_file.is_open()) {
  //   for (rbmd::Id i = 0; i < _RBE_P; ++i) {
  //     fac_file << i  << " "<<  h_fac[i]<<"\n";
  //   }
  //   fac_file.close();
  // }

  // 4. Calculate Rho for the *sampled* K-vectors (K_Sample)
  thrust::device_vector<rbmd::Real> d_rho_real(_RBE_P);
  thrust::device_vector<rbmd::Real> d_rho_imag(_RBE_P);
  thrust::device_vector<rbmd::Real> d_rho_npt_real(_RBE_P);
  thrust::device_vector<rbmd::Real> d_rho_npt_imag(_RBE_P);
  thrust::device_vector<rbmd::Real> d_energy_parts(2);  // [0]=sample, [1]=direct

  op::ComputePnumberChargeStructureFactorSOGOp<device::DEVICE_GPU>()(
      *_box, num_atoms, _RBE_P,
      thrust::raw_pointer_cast(_device_data->_d_charge.data()),
      thrust::raw_pointer_cast(_d_rbsog_K_Sample_x.data()),
      thrust::raw_pointer_cast(_d_rbsog_K_Sample_y.data()),
      thrust::raw_pointer_cast(_d_rbsog_K_Sample_z.data()),
      thrust::raw_pointer_cast(_device_data->_d_px.data()),
      thrust::raw_pointer_cast(_device_data->_d_py.data()),
      thrust::raw_pointer_cast(_device_data->_d_pz.data()),
      thrust::raw_pointer_cast(d_rho_real.data()),
      thrust::raw_pointer_cast(d_rho_imag.data()));

  // thrust::host_vector<rbmd::Real> h_rho_real = d_rho_real;
  // thrust::host_vector<rbmd::Real> h_rho_imag = d_rho_imag;
  // std::ofstream rho_file("h_rho.txt");
  // if (rho_file.is_open()) {
  //   for (rbmd::Id i = 0; i < _RBE_P; ++i) {
  //     rho_file << i  << " "<<  h_rho_real[i] << " "<< h_rho_imag[i]<<"\n";
  //   }
  //   rho_file.close();
  // }

  // 5. Calculate Force/Virial/Energy from *Sampled* K-vectors
  thrust::device_vector<rbmd::Real> d_SampleForce_x(num_atoms);
  thrust::device_vector<rbmd::Real> d_SampleForce_y(num_atoms);
  thrust::device_vector<rbmd::Real> d_SampleForce_z(num_atoms);
  thrust::device_vector<rbmd::Real> virial_sample(6);
  op::ComputeRBSOGSampleForceOp<device::DEVICE_GPU>()(
      *_box, num_atoms, _RBE_P, _qqr2e, _rbsog_S0, _rbsog_S_npt0,
      thrust::raw_pointer_cast(_d_rbsog_K_Sample_x.data()),
      thrust::raw_pointer_cast(_d_rbsog_K_Sample_y.data()),
      thrust::raw_pointer_cast(_d_rbsog_K_Sample_z.data()),
      thrust::raw_pointer_cast(_d_rbsog_K_npt_x.data()),
      thrust::raw_pointer_cast(_d_rbsog_K_npt_y.data()),
      thrust::raw_pointer_cast(_d_rbsog_K_npt_z.data()),
      thrust::raw_pointer_cast(_d_rbsog_idx_npt_all.data()),
      thrust::raw_pointer_cast(_d_fac.data()),
      thrust::raw_pointer_cast(_d_fac_npt.data()),
      thrust::raw_pointer_cast(d_rho_real.data()),
      thrust::raw_pointer_cast(d_rho_imag.data()),
      thrust::raw_pointer_cast(_device_data->_d_charge.data()),
      thrust::raw_pointer_cast(_device_data->_d_px.data()),
      thrust::raw_pointer_cast(_device_data->_d_py.data()),
      thrust::raw_pointer_cast(_device_data->_d_pz.data()),
      thrust::raw_pointer_cast(d_SampleForce_x.data()),
      thrust::raw_pointer_cast(d_SampleForce_y.data()),
      thrust::raw_pointer_cast(d_SampleForce_z.data()),
      thrust::raw_pointer_cast(virial_sample.data()),
      thrust::raw_pointer_cast(d_energy_parts.data()));

  // 6. Calculate Force/Virial/Energy from *Direct Sum* K-vectors
  // Host-side: generate the list of direct-sum K-vectors
  int bound = CEIL(SQRT(_rbsog_Kcut));
  rbmd::Real pxyz[3] = {(rbmd::Real)(2 * M_PI / xprd),
                        (rbmd::Real)(2 * M_PI / yprd),
                        (rbmd::Real)(2 * M_PI / zprd)};

  std::vector<rbmd::Real> h_k_direct_x, h_k_direct_y, h_k_direct_z;
  std::vector<rbmd::Real> h_f_b_sigma;
  std::vector<rbmd::Real> h_f_b_sigma_npt;

  for (int i = -bound; i <= bound; i++) {
    rbmd::Real Kx = pxyz[0] * (i + 0.00);
    for (int j = -bound; j <= bound; j++) {
      rbmd::Real Ky = pxyz[1] * (j + 0.00);
      for (int k = -bound; k <= bound; k++) {
        rbmd::Real Kz = pxyz[2] * (k + 0.00);
        if ((!(i == 0 && j == 0 && k == 0)) &&
            (i * i + j * j + k * k <= _rbsog_Kcut)) {
          h_k_direct_x.push_back(Kx);
          h_k_direct_y.push_back(Ky);
          h_k_direct_z.push_back(Kz);
          h_f_b_sigma.push_back(
              Gaussian_Fourier_Plus(Kx, Ky, Kz, _rbsog_sigma, _rbsog_b,
                                    _rbsog_w0, _rbsog_Mmax, _h_rbsog_coef));

          h_f_b_sigma_npt.push_back(Gaussian_Fourier_Plus_modify(
              Kx, Ky, Kz, _rbsog_sigma, _rbsog_b, _rbsog_w0, _rbsog_Mmax,
              _h_rbsog_coef_npt));
        }
      }
    }
  }



  auto num_k_direct = h_k_direct_x.size();

  // for (int i = 0; i < num_k_direct; ++i) {
  //   std::cout << "h_k_direct: "<< h_k_direct_x[i] <<
  //     ", "<<h_k_direct_y[i] << ", " <<  h_k_direct_z[i]  <<std::endl;
  // }

  // Copy direct-sum data to device
  thrust::device_vector<rbmd::Real> d_k_direct_x = h_k_direct_x;
  thrust::device_vector<rbmd::Real> d_k_direct_y = h_k_direct_y;
  thrust::device_vector<rbmd::Real> d_k_direct_z = h_k_direct_z;

  thrust::device_vector<rbmd::Real> d_f_b_sigma = h_f_b_sigma;
  thrust::device_vector<rbmd::Real> d_f_b_sigma_npt = h_f_b_sigma_npt;

  // Calculate Rho for direct-sum vectors
  thrust::device_vector<rbmd::Real> d_rho_direct_real(num_k_direct);
  thrust::device_vector<rbmd::Real> d_rho_direct_imag(num_k_direct);
  op::ComputeDirectChargeStructureFactorOp<device::DEVICE_GPU>()(
      num_atoms, num_k_direct, thrust::raw_pointer_cast(d_k_direct_x.data()),
      thrust::raw_pointer_cast(d_k_direct_y.data()),
      thrust::raw_pointer_cast(d_k_direct_z.data()),
      thrust::raw_pointer_cast(_device_data->_d_charge.data()),
      thrust::raw_pointer_cast(_device_data->_d_px.data()),
      thrust::raw_pointer_cast(_device_data->_d_py.data()),
      thrust::raw_pointer_cast(_device_data->_d_pz.data()),
      thrust::raw_pointer_cast(d_rho_direct_real.data()),
      thrust::raw_pointer_cast(d_rho_direct_imag.data()));

  // Calculate Force from direct-sum vectors
  thrust::device_vector<rbmd::Real> d_DirectForce_x(num_atoms);
  thrust::device_vector<rbmd::Real> d_DirectForce_y(num_atoms);
  thrust::device_vector<rbmd::Real> d_DirectForce_z(num_atoms);
  thrust::device_vector<rbmd::Real> virial_direct(6);
  op::ComputeRBSOGDirectForceOp<device::DEVICE_GPU>()(
      *_box, num_atoms, num_k_direct, _qqr2e,
      thrust::raw_pointer_cast(d_k_direct_x.data()),
      thrust::raw_pointer_cast(d_k_direct_y.data()),
      thrust::raw_pointer_cast(d_k_direct_z.data()),
      thrust::raw_pointer_cast(d_f_b_sigma.data()),
      thrust::raw_pointer_cast(d_f_b_sigma_npt.data()),
      thrust::raw_pointer_cast(d_rho_direct_real.data()),
      thrust::raw_pointer_cast(d_rho_direct_imag.data()),
      thrust::raw_pointer_cast(_device_data->_d_charge.data()),
      thrust::raw_pointer_cast(_device_data->_d_px.data()),
      thrust::raw_pointer_cast(_device_data->_d_py.data()),
      thrust::raw_pointer_cast(_device_data->_d_pz.data()),
      thrust::raw_pointer_cast(d_DirectForce_x.data()),
      thrust::raw_pointer_cast(d_DirectForce_y.data()),
      thrust::raw_pointer_cast(d_DirectForce_z.data()),
      thrust::raw_pointer_cast(virial_direct.data()),
      thrust::raw_pointer_cast(d_energy_parts.data()));


  TransformForces(_device_data->_d_force_kspace_x, d_SampleForce_x,
                  d_DirectForce_x);
  TransformForces(_device_data->_d_force_kspace_y, d_SampleForce_y,
                  d_DirectForce_y);
  TransformForces(_device_data->_d_force_kspace_z, d_SampleForce_z,
                  d_DirectForce_z);

  //virial      Virial_Total = Virial_sample + Virial_direct
  TransformForces(_device_data->_d_virial_kspace, virial_sample,
                virial_direct);

  // thrust::host_vector<rbmd::Real> h_kspace_virial =_device_data->_d_virial_kspace;
  //
  // std::ofstream kspace_file("kspace_rbsog_virial.txt");
  // if (kspace_file.is_open()) {
  //   for (rbmd::Id i = 0; i < h_kspace_virial.size(); ++i) {
  //     kspace_file << i  << " " <<h_kspace_virial[i]  << "\n";
  //   }
  //   kspace_file.close();
  // }

  // thrust::host_vector<rbmd::Real> h_kspace_x =_device_data->_d_force_kspace_x;
  // thrust::host_vector<rbmd::Real> h_kspace_y = _device_data->_d_force_kspace_y;
  // thrust::host_vector<rbmd::Real> h_kspace_z = _device_data->_d_force_kspace_z;
  //
  // std::ofstream kspace_file("kspace_rbsog.txt");
  // auto atom_id_to_idx =
  //   LinkedCellLocator::GetInstance().GetLinkedCell()->_atom_id_to_idx;
  // if (kspace_file.is_open()) {
  //   for (rbmd::Id i = 0; i < h_kspace_x.size(); ++i) {
  //     auto index = atom_id_to_idx[i];
  //     kspace_file << i  << " " <<h_kspace_x[index] <<" " <<h_kspace_y[index] <<" " <<
  //     h_kspace_z[index]<< "\n";
  //   }
  //   kspace_file.close();
  // }

  // // 9. Finalize Energy
   thrust::host_vector<rbmd::Real> h_energy_parts = d_energy_parts;
   rbmd::Real  energy_sample = h_energy_parts[0];
   rbmd::Real  energy_direct = h_energy_parts[1];

  rbmd::Real coeff = (LOG(_rbsog_b) / (SQRT(2 * M_PI) * _rbsog_sigma)) *
    (_rbsog_w0 + (1 - POW(_rbsog_b, -_rbsog_Mmax)) / (_rbsog_b - 1));

  rbmd::Real self_energy = -coeff * _sum_sq_charge * _qqr2e;
  _e_kspace = energy_sample + energy_direct + self_energy;


  // 10. Timing and Virial Reduction
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration = end - start;
  TimingStatistics::Instance().record("Long-Range", duration.count());

}

rbmd::Real KSpaceCalculator::G_sigma(rbmd::Real sigma, rbmd::Real r) const {
  return EXP(-r * r / (2 * sigma * sigma)) / SQRT(2 * M_PI * sigma * sigma);//标准高斯函数
}

rbmd::Real KSpaceCalculator::Compute_W0(rbmd::Real r0, rbmd::Real b) const {
  rbmd::Real sum = 0.00;
  for (int i = 1; i < 200; i++) {
    sum = sum + POW(b, (rbmd::Real)(-i)) *
                    G_sigma(1.0, POW(b, (rbmd::Real)(-i)) * r0);
  }
  rbmd::Real w0 = (1.0 / G_sigma(1.0, r0)) * ((1.0 / (2 * LOG(b) * r0)) - sum);
  return w0;
}

rbmd::Real KSpaceCalculator::Gaussian(
    int kx, int ky, int kz, const Box& box, rbmd::Real sigma, rbmd::Real b,
    rbmd::Real w0, int Mmax,
    const thrust::host_vector<rbmd::Real>& coef) const {
  rbmd::Real k2 = POW((kx + 0.00) * 2 * M_PI / box._length[0], 2) +
                  POW((ky + 0.00) * 2 * M_PI / box._length[1], 2) +
                  POW((kz + 0.00) * 2 * M_PI / box._length[2], 2);
  rbmd::Real sum = 0.00;
  rbmd::Real b2 = b * b;
  for (int i = 0; i < Mmax; i++) {
    sum = sum + coef[i] * EXP(-(POW(b2, (i + 0.0)) * sigma * sigma) * k2 / 2.0);
  }
  return sum;
}

rbmd::Real KSpaceCalculator::Gaussian_modify(
    int kx, int ky, int kz, const Box& box, rbmd::Real sigma, rbmd::Real b,
    rbmd::Real w0, int Mmax,
    const thrust::host_vector<rbmd::Real>& coef_npt) const {
  rbmd::Real k2 = POW((kx + 0.00) * 2 * M_PI / box._length[0], 2) +
                  POW((ky + 0.00) * 2 * M_PI / box._length[1], 2) +
                  POW((kz + 0.00) * 2 * M_PI / box._length[2], 2);
  rbmd::Real sum = 0.00;
  rbmd::Real b2 = b * b;
  for (int i = 0; i < Mmax; i++) {
    sum = sum +
          coef_npt[i] * EXP(-(POW(b2, (i + 0.0)) * sigma * sigma) * k2 / 2.0);
  }
  return sum;
}

rbmd::Real KSpaceCalculator::Gaussian_Fourier_Plus(
    rbmd::Real Kx, rbmd::Real Ky, rbmd::Real Kz, rbmd::Real sigma, rbmd::Real b,
    rbmd::Real w0, rbmd::Id Mmax,
    const thrust::host_vector<rbmd::Real>& coef) const {
  rbmd::Real k2 = Kx * Kx + Ky * Ky + Kz * Kz;
  rbmd::Real sum = 0.00;
  rbmd::Real b2 = b * b;
  for (int i = 0; i < Mmax; i++) {
    sum = sum + coef[i] * EXP(-(POW(b2, (i + 0.0)) * sigma * sigma) * k2 / 2.0);
  }
  return sum;
}

rbmd::Real KSpaceCalculator::Gaussian_Fourier_Plus_modify(
    rbmd::Real Kx, rbmd::Real Ky, rbmd::Real Kz, rbmd::Real sigma, rbmd::Real b,
    rbmd::Real w0, int Mmax,
    const thrust::host_vector<rbmd::Real>& coef_npt) const {
  rbmd::Real k2 = Kx * Kx + Ky * Ky + Kz * Kz;
  rbmd::Real sum = 0.00;
  rbmd::Real b2 = b * b;
  for (int i = 0; i < Mmax; i++) {
    sum = sum +
          coef_npt[i] * EXP(-(POW(b2, (i + 0.0)) * sigma * sigma) * k2 / 2.0);
  }
  return sum;
}

rbmd::Real KSpaceCalculator::randn_box_muller(rbmd::Real Mean,
                                              rbmd::Real SquareMargin) {
  const rbmd::Real epsilon = 1.17549e-038;
  const rbmd::Real two_pi = 2.0 * M_PI;

  // std::uniform_real_distribution<rbmd::Real> dis(0.0, 1.0);
  rbmd::Real u1, u2;
  do {
    u1 = RandomValue<rbmd::Real>(0.0, 1.0);
    u2 = RandomValue<rbmd::Real>(0.0, 1.0);
  } while (u1 <= epsilon);

  rbmd::Real z0 = SQRT(-2.0 * LOG(u1)) * COS(two_pi * u2);
  rbmd::Real z1 = SQRT(-2.0 * LOG(u1)) * SIN(two_pi * u2);

  return z0 * SquareMargin + Mean;
}

rbmd::Real KSpaceCalculator::MH_D_Modify(int xx, rbmd::Real factor) const {
  rbmd::Real qes;
  if (xx == 0)
    qes = ERF(0.5 * factor);
  else
    qes = 0.5 * (ERF((ABS(xx) + 0.5) * factor) - ERF((ABS(xx) - 0.5) * factor));

  return qes;
}
