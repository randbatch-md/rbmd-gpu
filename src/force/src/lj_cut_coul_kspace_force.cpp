#include "lj_cut_coul_kspace_force.h"

#include <thrust/device_ptr.h>

#include "../../common/device_types.h"
#include "../../common/rbmd_define.h"
#include "../../common/types.h"
#include "ljforce_op/ljforce_op.h"
#include "../common/RBEPSample.h"
#include "neighbor_list/include/neighbor_list_builder/half_neighbor_list_builder.h"
#include "neighbor_list/include/neighbor_list_builder/full_neighbor_list_builder.h"
#include "neighbor_list/include/neighbor_list_builder/rbl_full_neighbor_list_builder.h"
// #include <hipcub/hipcub.hpp>
// #include <hipcub/backend/rocprim/block/block_reduce.hpp>

extern int test_current_step;

LJCutCoulKspaceForce::LJCutCoulKspaceForce()
{
  _rbl_neighbor_list_builder = std::make_shared<RblFullNeighborListBuilder>();
  _neighbor_list_builder = std::make_shared<FullNeighborListBuilder>();

  CHECK_RUNTIME(MALLOC(&_d_total_evdwl, sizeof(rbmd::Real)));
  CHECK_RUNTIME(MALLOC(&_d_total_ecoul, sizeof(rbmd::Real)));
}

LJCutCoulKspaceForce::~LJCutCoulKspaceForce()
{
  CHECK_RUNTIME(FREE(_d_total_evdwl));
  CHECK_RUNTIME(FREE(_d_total_ecoul));

}

void LJCutCoulKspaceForce::Init()
{
    _num_atoms = *(_structure_info_data->_num_atoms);
    _corr_value_x = 0;
    _corr_value_y = 0;
    _corr_value_z = 0;

     _ave_evdwl = 0.0;
     _ave_ecoul = 0.0;
     _ave_self_energy = 0.0;
     _ave_ekspace = 0.0;
     _qqr2e = 1.0;

   _cut_off = DataManager::getInstance().getConfigData()->Get
    <rbmd::Real>("cut_off", "hyper_parameters", "neighbor");
   _neighbor_type = DataManager::getInstance().getConfigData()->Get
      <std::string>("type", "hyper_parameters", "neighbor");

   _coulomb_type =DataManager::getInstance().getConfigData()->Get<std::string>(
        "type", "hyper_parameters", "coulomb");
    _RBE_P = DataManager::getInstance().getConfigData()->Get<rbmd::Id>(
      "coulomb_sample_num", "hyper_parameters", "coulomb");
    _alpha = DataManager::getInstance().getConfigData()->Get<rbmd::Real>(
  "alpha", "hyper_parameters", "coulomb");
    _Kmax = DataManager::getInstance().getConfigData()->Get<rbmd::Id>(
"kmax", "hyper_parameters", "coulomb");
    _num_k =  POW(2 * _Kmax + 1,3.0) - 1;

    RBEInit(_device_data->_d_box,_alpha,_RBE_P);
}

void LJCutCoulKspaceForce::Execute()
{
  ComputeLJCutCoulForce();
  ComputeKspaceForce();
  SumForces();
}

void LJCutCoulKspaceForce::ComputeLJCutCoulForce()
{

  //
  if ("RBL" ==_neighbor_type)
  {
    // rbl_neighbor_list_build
    auto start = std::chrono::high_resolution_clock::now();
    rbl_list = _rbl_neighbor_list_builder->Build();

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

    op::LJCutCoulRBLForceOp<device::DEVICE_GPU> lj_cut_coul_rbl_force_op;
    lj_cut_coul_rbl_force_op(
        _device_data->_d_box, r_core, _cut_off, _num_atoms,
        neighbor_sample_num,rbl_list->_selection_frequency,_alpha,_qqr2e,
        thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
        thrust::raw_pointer_cast(_device_data->_d_molecular_id.data()),
        thrust::raw_pointer_cast(_device_data->_d_sigma.data()),
        thrust::raw_pointer_cast(_device_data->_d_eps.data()),
        thrust::raw_pointer_cast(rbl_list->_start_idx.data()),
        thrust::raw_pointer_cast(rbl_list->_end_idx.data()),
        thrust::raw_pointer_cast(rbl_list->_d_neighbors.data()),
        thrust::raw_pointer_cast(rbl_list->_d_random_neighbor.data()),
        thrust::raw_pointer_cast(rbl_list->_d_random_neighbor_num.data()),
        thrust::raw_pointer_cast(_device_data->_d_charge.data()),
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

    // energy
    list = _neighbor_list_builder->Build();

    rbmd::Real h_total_evdwl = 0.0;
    rbmd::Real h_total_ecoul = 0.0;

    CHECK_RUNTIME(MEMSET(_d_total_evdwl, 0, sizeof(rbmd::Real)));
    CHECK_RUNTIME(MEMSET(_d_total_ecoul, 0, sizeof(rbmd::Real)));

    op::LJCutCoulEnergyOp<device::DEVICE_GPU> lj_cut_coul_energy_op;
    lj_cut_coul_energy_op(_device_data->_d_box, _cut_off, _num_atoms,_alpha,_qqr2e,
                  thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
                  thrust::raw_pointer_cast(_device_data->_d_molecular_id.data()),
                  thrust::raw_pointer_cast(_device_data->_d_sigma.data()),
                  thrust::raw_pointer_cast(_device_data->_d_eps.data()),
                  thrust::raw_pointer_cast(list->_start_idx.data()),
                  thrust::raw_pointer_cast(list->_end_idx.data()),
                  thrust::raw_pointer_cast(list->_d_neighbors.data()),
                  thrust::raw_pointer_cast(_device_data->_d_charge.data()),
                  thrust::raw_pointer_cast(_device_data->_d_px.data()),
                  thrust::raw_pointer_cast(_device_data->_d_py.data()),
                  thrust::raw_pointer_cast(_device_data->_d_pz.data()),
                  _d_total_evdwl,_d_total_ecoul);

    CHECK_RUNTIME(MEMCPY(&h_total_evdwl,_d_total_evdwl , sizeof(rbmd::Real), D2H));
    CHECK_RUNTIME(MEMCPY(&h_total_ecoul,_d_total_ecoul , sizeof(rbmd::Real), D2H));

    // 打印累加后的总能量
    _ave_evdwl = h_total_evdwl/_num_atoms;
    _ave_ecoul = h_total_ecoul/_num_atoms;

    std::cout << "test_current_step:" << test_current_step <<  " ,"
    << "average_vdwl_energy:" << _ave_evdwl << " ," <<  "average_coul_energy:" << _ave_ecoul << std::endl;

    //out
    std::ofstream outfile("ave_ljcoul_rbl.txt", std::ios::app);
    outfile << test_current_step << " " << _ave_evdwl  << " "<< _ave_ecoul << std::endl;
    outfile.close();
  }
  else
  {
      //neighbor_list_build
  auto start = std::chrono::high_resolution_clock::now();
  list = _neighbor_list_builder->Build();

  auto end = std::chrono::high_resolution_clock::now();

  std::chrono::duration<rbmd::Real> duration = end - start;
  std::cout << "构建verlet-list耗时" << duration.count() << "秒" << std::endl;

  //

  rbmd::Real h_total_evdwl = 0.0;
  rbmd::Real h_total_ecoul = 0.0;

  CHECK_RUNTIME(MEMSET(_d_total_evdwl, 0, sizeof(rbmd::Real)));
  CHECK_RUNTIME(MEMSET(_d_total_ecoul, 0, sizeof(rbmd::Real)));

  //
  op::LJCutCoulForceOp<device::DEVICE_GPU> lj_cut_coul_force_op;
  lj_cut_coul_force_op(_device_data->_d_box, _cut_off, _num_atoms,_alpha,_qqr2e,
                  thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
                  thrust::raw_pointer_cast(_device_data->_d_molecular_id.data()),
                  thrust::raw_pointer_cast(_device_data->_d_sigma.data()),
                  thrust::raw_pointer_cast(_device_data->_d_eps.data()),
                  thrust::raw_pointer_cast(list->_start_idx.data()),
                  thrust::raw_pointer_cast(list->_end_idx.data()),
                  thrust::raw_pointer_cast(list->_d_neighbors.data()),
                  thrust::raw_pointer_cast(_device_data->_d_charge.data()),
                  thrust::raw_pointer_cast(_device_data->_d_px.data()),
                  thrust::raw_pointer_cast(_device_data->_d_py.data()),
                  thrust::raw_pointer_cast(_device_data->_d_pz.data()),
                  thrust::raw_pointer_cast(_device_data->_d_force_ljcoul_x.data()),
                  thrust::raw_pointer_cast(_device_data->_d_force_ljcoul_y.data()),
                  thrust::raw_pointer_cast(_device_data->_d_force_ljcoul_z.data()),
                  _d_total_evdwl,_d_total_ecoul);

  CHECK_RUNTIME(MEMCPY(&h_total_evdwl,_d_total_evdwl , sizeof(rbmd::Real), D2H));
  CHECK_RUNTIME(MEMCPY(&h_total_ecoul,_d_total_ecoul , sizeof(rbmd::Real), D2H));

  // 打印累加后的总能量
  _ave_evdwl = h_total_evdwl/_num_atoms;
  _ave_ecoul = h_total_ecoul/_num_atoms;

  std::cout << "test_current_step:" << test_current_step <<  " ,"
  << "average_vdwl_energy:" << _ave_evdwl << " ," <<  "average_coul_energy:" << _ave_ecoul << std::endl;

  //out
  std::ofstream outfile("ave_ljcoul.txt", std::ios::app);
  outfile << test_current_step << " " << _ave_evdwl  << " "<< _ave_ecoul << std::endl;
  outfile.close();
  }
}

void LJCutCoulKspaceForce::ComputeKspaceForce()
{
  if("RBE" == _coulomb_type)
  {
      ComputeRBEForce();
  }
  else
  {
     ComputeEwladForce();
  }
}

void LJCutCoulKspaceForce::SumForces()
{

  thrust::transform(
      _device_data->_d_force_ljcoul_x.begin(), _device_data->_d_force_ljcoul_x.end(),
      _device_data->_d_force_ewald_x.begin(),
      _device_data->_d_fx.begin(), // 将结果存回到 _d_fx 中
      thrust::plus<rbmd::Real>());

  thrust::transform(
      _device_data->_d_force_ljcoul_y.begin(), _device_data->_d_force_ljcoul_y.end(),
      _device_data->_d_force_ewald_y.begin(),
      _device_data->_d_fy.begin(), // 将结果存回到 _d_fy 中
      thrust::plus<rbmd::Real>());

  thrust::transform(
      _device_data->_d_force_ljcoul_z.begin(), _device_data->_d_force_ljcoul_z.end(),
      _device_data->_d_force_ewald_z.begin(),
      _device_data->_d_fz.begin(), // 将结果存回到 _d_fz 中
      thrust::plus<rbmd::Real>());
}

void LJCutCoulKspaceForce::ComputeChargeStructureFactorEwald(
    Box* box,
    rbmd::Id num_atoms,
    rbmd::Id Kmax,
    rbmd::Real alpha,
    rbmd::Real qqr2e,
    rbmd::Real* value_Re_array,
    rbmd::Real* value_Im_array)
{
    //thrust::fill(density_real.begin(), density_real.end(), 0.0f);
    //thrust::fill(density_imag.begin(), density_imag.end(), 0.0f);
    thrust::device_vector<rbmd::Real> density_real_atom;
    thrust::device_vector<rbmd::Real> density_imag_atom;
    density_real_atom.resize(_num_atoms);
    density_imag_atom.resize(_num_atoms);

    rbmd::Real total_energy_ewald = 0;
    rbmd::Id index = 0;
    for (rbmd::Id i = -Kmax; i <= Kmax; i++)
    {
        for (rbmd::Id j = -Kmax; j <= Kmax; j++)
        {
            for (rbmd::Id k = -Kmax; k <= Kmax; k++)
            {
                if (!(i == 0 && j == 0 && k == 0))
                {
                    Real3 K = make_Real3(rbmd::Real(2 * M_PI * i / box->_length[0]),
                                           rbmd::Real(2 * M_PI * j / box->_length[1]),
                                           rbmd::Real(2 * M_PI * k / box->_length[2]));
                    rbmd::Real Range_K = SQRT(K.x * K.x + K.y * K.y + K.z * K.z);
                    rbmd::Real Range_K2 = Range_K*Range_K;
                    rbmd::Real alpha_inv =  1 / alpha;

                    op::ComputeChargeStructureFactorOp<device::DEVICE_GPU> charge_structure_factor_op;
                    charge_structure_factor_op(num_atoms, K,
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

                    value_Re_array[index] = value_Re;
                    value_Im_array[index] = value_Im;
                    index++;
                }
            }
        }
    }

  //charge self energy//
  ComputeSelfEnergy(alpha,qqr2e,_ave_self_energy);

  //compute ewald energy//
  rbmd::Real volume = box->_length[0] * box->_length[1]*box->_length[2];
  total_energy_ewald = qqr2e * (2 * M_PI / volume) * total_energy_ewald;
  _ave_ekspace = total_energy_ewald / _num_atoms;

  _ave_ekspace = _ave_ekspace + _ave_self_energy;

  //out
   std::cout << "test_current_step:" << test_current_step <<  " ,"
   << "ave_eewald:" << _ave_ekspace << std::endl;

  std::ofstream outfile("ave_ewald.txt", std::ios::app);
  outfile << test_current_step << " "<< _ave_ekspace << std::endl;
  outfile.close();
}

void LJCutCoulKspaceForce::ComputeEwladForce()
{
    rbmd::Real* value_Re_array;
    rbmd::Real* value_Im_array;

    //EwaldForce//
    CHECK_RUNTIME(MALLOC(&value_Re_array, _num_k * sizeof(rbmd::Real)));
    CHECK_RUNTIME(MALLOC(&value_Im_array, _num_k * sizeof(rbmd::Real)));
    CHECK_RUNTIME(MEMSET(value_Re_array, 0, _num_k *sizeof(rbmd::Real)));
    CHECK_RUNTIME(MEMSET(value_Im_array, 0, _num_k *sizeof(rbmd::Real)));

    ComputeChargeStructureFactorEwald(_device_data->_d_box, _num_atoms, _Kmax,
      _alpha,_qqr2e, value_Re_array, value_Im_array);


    op::ComputeEwaldForceOp<device::DEVICE_GPU> ewlad_force_op;
    ewlad_force_op(_device_data->_d_box,_num_atoms, _Kmax, _alpha,_qqr2e,
        value_Re_array,value_Im_array,
        thrust::raw_pointer_cast(_device_data->_d_charge.data()),
        thrust::raw_pointer_cast(_device_data->_d_px.data()),
        thrust::raw_pointer_cast(_device_data->_d_py.data()),
        thrust::raw_pointer_cast(_device_data->_d_pz.data()),
        thrust::raw_pointer_cast(_device_data->_d_force_ewald_x.data()),
        thrust::raw_pointer_cast(_device_data->_d_force_ewald_y.data()),
        thrust::raw_pointer_cast(_device_data->_d_force_ewald_z.data()));

    CHECK_RUNTIME(FREE(value_Re_array));
    CHECK_RUNTIME(FREE(value_Im_array));
}

void LJCutCoulKspaceForce::RBEInit(Box* box,rbmd::Real alpha,rbmd::Id RBE_P)
{
  Real3 sigma = { rbmd::Real((SQRT(alpha / 2.0) * box->_length[0]/M_PI)),
                  rbmd::Real((SQRT(alpha / 2.0) * box->_length[1]/M_PI)),
                  rbmd::Real((SQRT(alpha / 2.0) * box->_length[2]/M_PI))};
  auto random = true;
  RBEPSAMPLE rbe_presolve_psample = { alpha, box, RBE_P, random};

  _P_Sample_x.resize(RBE_P);
  _P_Sample_y.resize(RBE_P);
  _P_Sample_z.resize(RBE_P);

  rbe_presolve_psample.Fetch_P_Sample(0.0, sigma,
    thrust::raw_pointer_cast(_P_Sample_x.data()),
    raw_pointer_cast(_P_Sample_y.data()),
    raw_pointer_cast(_P_Sample_z.data()));

  //index key
  _psample_key.resize(_num_atoms * RBE_P);
  op::GenerateIndexArrayOp<device::DEVICE_GPU> generate_index_array_op;
  generate_index_array_op(_num_atoms,RBE_P,
    thrust::raw_pointer_cast(_psample_key.data()));
}

void LJCutCoulKspaceForce::ComputeChargeStructureFactorRBE(
   Box* box,
   rbmd::Id num_atoms,
   rbmd::Id Kmax,
   rbmd::Real alpha,
   rbmd::Id RBE_P,
   rbmd::Real qqr2e,
   thrust::device_vector<rbmd::Real> rhok_real_redue,
   thrust::device_vector<rbmd::Real> rhok_image_redue)
{
  //
  thrust::device_vector<rbmd::Real>  rhok_real_atom;
  thrust::device_vector<rbmd::Real>  rhok_image_atom;

  rhok_real_atom.resize(_num_atoms* RBE_P);
  rhok_image_atom.resize(_num_atoms* RBE_P);
  auto p_number= RBE_P;

  //Charge Structure Factor
  op::ComputePnumberChargeStructureFactorOp<device::DEVICE_GPU> pnumber_charge_structure_factor_op;
  pnumber_charge_structure_factor_op(_device_data->_d_box, _num_atoms, p_number,
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


  //charge self energy//
  ComputeSelfEnergy(alpha,qqr2e,_ave_self_energy);

  //ave_eewald
  ComputeEwaldEnergy(_device_data->_d_box, _num_atoms, Kmax,
      alpha, qqr2e ,_ave_ekspace);
  _ave_ekspace = _ave_ekspace +_ave_self_energy;

    //out
   std::cout << "test_current_step:" << test_current_step <<  " ,"
   << "ave_ekspace:" << _ave_ekspace << std::endl;

  std::ofstream outfile("ave_rbe.txt", std::ios::app);
  outfile << test_current_step << " "<< _ave_ekspace << std::endl;
  outfile.close();
}

void LJCutCoulKspaceForce::ComputeRBEForce()
{
  //
  _rhok_real_redue.resize(_RBE_P);
  _rhok_image_redue.resize(_RBE_P);
  ComputeChargeStructureFactorRBE(_device_data->_d_box, _num_atoms, _Kmax,
      _alpha,_RBE_P,_qqr2e,_rhok_real_redue,_rhok_image_redue);

   //RBE Force
  op::ComputeRBEForceOp<device::DEVICE_GPU> rbe_force_op;
  rbe_force_op(_device_data->_d_box,_num_atoms, _RBE_P,_alpha,_qqr2e,
        thrust::raw_pointer_cast(_rhok_real_redue.data()),
        thrust::raw_pointer_cast(_rhok_image_redue.data()),
        thrust::raw_pointer_cast(_device_data->_d_charge.data()),
        thrust::raw_pointer_cast(_P_Sample_x.data()),
        raw_pointer_cast(_P_Sample_y.data()),
        raw_pointer_cast(_P_Sample_z.data()),
        thrust::raw_pointer_cast(_device_data->_d_px.data()),
        thrust::raw_pointer_cast(_device_data->_d_py.data()),
        thrust::raw_pointer_cast(_device_data->_d_pz.data()),
        thrust::raw_pointer_cast(_device_data->_d_force_ewald_x.data()),
        thrust::raw_pointer_cast(_device_data->_d_force_ewald_y.data()),
        thrust::raw_pointer_cast(_device_data->_d_force_ewald_z.data()));
}

void LJCutCoulKspaceForce::ComputeSelfEnergy(
  rbmd::Real  alpha,
  rbmd::Real  qqr2e,
  rbmd::Real& ave_self_energy)
{
  //compute self_energy
  thrust::device_vector<rbmd::Real> sq_charge;
  sq_charge.resize(_num_atoms);

  op::SqchargeOp<device::DEVICE_GPU> sq_charge_op;
  sq_charge_op(_num_atoms,
    thrust::raw_pointer_cast(_device_data->_d_charge.data()),
    thrust::raw_pointer_cast(sq_charge.data()));

  rbmd::Real sum_sq_charge = thrust::reduce(sq_charge.begin(), sq_charge.end(), 0.0f, thrust::plus<rbmd::Real>());
  rbmd::Real total_self_energy = qqr2e * (- sqrt(alpha / M_PI) *sum_sq_charge);

  ave_self_energy =  total_self_energy / _num_atoms;
}

void LJCutCoulKspaceForce::ComputeEwaldEnergy(
    Box* box,
    rbmd::Id num_atoms,
    rbmd::Id Kmax,
    rbmd::Real alpha,
    rbmd::Real qqr2e,
    rbmd::Real&  ave_ekspace)
{
    thrust::device_vector<rbmd::Real> density_real_atom;
    thrust::device_vector<rbmd::Real> density_imag_atom;
    density_real_atom.resize(num_atoms);
    density_imag_atom.resize(num_atoms);

    rbmd::Real total_energy_ewald = 0;
    for (rbmd::Id i = -Kmax; i <= Kmax; i++)
    {
        for (rbmd::Id j = -Kmax; j <= Kmax; j++)
        {
            for (rbmd::Id k = -Kmax; k <= Kmax; k++)
            {
                if (!(i == 0 && j == 0 && k == 0))
                {
                    Real3 K = make_Real3(rbmd::Real(2 * M_PI * i / box->_length[0]),
                                           rbmd::Real(2 * M_PI * j / box->_length[1]),
                                           rbmd::Real(2 * M_PI * k / box->_length[2]));
                    rbmd::Real Range_K = SQRT(K.x * K.x + K.y * K.y + K.z * K.z);
                    rbmd::Real Range_K2 = Range_K*Range_K;
                    rbmd::Real alpha_inv =  1 / alpha;

                    op::ComputeChargeStructureFactorOp<device::DEVICE_GPU> charge_structure_factor_op;
                    charge_structure_factor_op(num_atoms, K,
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

  rbmd::Real volume = box->_length[0] * box->_length[1]*box->_length[2];
  total_energy_ewald = qqr2e * (2 * M_PI / volume) * total_energy_ewald;
  ave_ekspace = total_energy_ewald / num_atoms;
}


