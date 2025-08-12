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
#include "output/include/Logger.hpp"

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

void KSpaceCalculator::Execute()
{
    if ("RBE" == _coulomb_type)
    {
        ComputeRBE();
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

    //
    ReduceVirial(num_atoms, _device_data->_d_flat_virial_kspace,
      _device_data->_d_virial_kspace);
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

  ComputeChargeStructureFactorRBE(*_box, num_atoms, _kmax_array,
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

  ComputeRBEVirial();

  //   ReduceVirial(num_atoms,_device_data->_d_flat_virial_kspace,
  // _device_data->_d_virial_kspace);
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

  _sum_charge = thrust::reduce(_device_data->_d_charge.begin(),
_device_data->_d_charge.end(), 0.0f, thrust::plus<rbmd::Real>());
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