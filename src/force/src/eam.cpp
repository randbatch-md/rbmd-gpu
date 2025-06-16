#include "eam.h"

#include <output/include/Logger.hpp>

#include "../../common/device_types.h"
#include "../../common/rbmd_define.h"
#include "../../common/types.h"
#include "../common/unit_factor.h"
#include "eam_op/eam_op.h"
#include "force_op/force_op.h"
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

EAM::EAM() {
  _rbl_neighbor_list_builder = std::make_shared<RblFullNeighborListBuilder>();
  _neighbor_list_builder = std::make_shared<FullNeighborListBuilder>();

  auto unit = DataManager::getInstance().getConfigData()->Get
  <std::string>("unit", "init_configuration", "read_data");
  UNIT unit_factor = unit_factor_map[unit];
  std::remove("thermo.txt");
}

EAM::~EAM() {  }

void EAM::Init() {
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
      Logger::Instance().info( "\033[31mFATAL ERROR: When using RBL for the neighbor type, "
                   "the key 'energy_rbl_flag' must be defined.\033[0m");
      exit(EXIT_FAILURE); //
    }
  }

  //
  std::string potential_file =
      DataManager::getInstance().getConfigData()->Get<std::string>(
          "potential_file", "hyper_parameters", "force_field");
  ReadPotentialFile(potential_file);
  InitStyle();
}

void EAM::Execute() {
  SumForces();
  EvaluatePotentialenergy();
}

void EAM::ReadPotentialFile(const std::string& filename) {
  std::ifstream input_file(filename);
  if (!input_file.is_open())
  {
    std::cerr << "Unable to open the file." << std::endl;
  }
  for (int i = 0; i < 2; ++i)
  {
    input_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
  input_file >> file.nrho >> file.drho >> file.nr >> file.dr >> file.cut_off;

  file.frho.resize(file.nrho + 1);
  file.zr.resize(file.nr + 1);
  file.rhor.resize(file.nrho + 1);

  for (int i = 0; i < file.nrho; ++i)
  {
    input_file >> file.frho[i];
  }

  for (int i = 0; i < file.nr; ++i)
  {
    input_file >> file.zr[i];
  }

  for (int i = 0; i < file.nrho; ++i)
  {
    input_file >> file.rhor[i];
  }
  std::cout<< "test--EAM--------" << file.nrho  << ", "<< file.drho << ", " << ", "
    <<file.nr   << ", " << file.dr << ", "  << file.cut_off  <<std::endl;

  input_file.close();
}

void EAM::InitStyle() {
  AllocateEAM();
  file2array();
  array2spline();
  SetEAM();
}

void EAM::AllocateEAM() {  }

void EAM::file2array() {
  rbmd::Id i, j, k, m, n;
  rbmd::Real sixth = 1.0 / 6.0;

  rbmd::Real rmax;
  eam_paras.dr = eam_paras.drho = rmax = eam_paras.rhomax = 0.0;

  eam_paras.dr = MAX(eam_paras.dr, file.dr);
  eam_paras.drho = MAX(eam_paras.drho, file.drho);
  rmax = MAX(rmax, (file.nr - 1) * file.dr);
  eam_paras.rhomax = MAX(eam_paras.rhomax, (file.nrho - 1) * file.drho);

  // set nr,nrho from cutoff and spacings
  // 0.5 is for round-off in divide

  eam_paras.nr = static_cast<int>(rmax / eam_paras.dr + 0.5);
  eam_paras.nrho = static_cast<int>(eam_paras.rhomax / eam_paras.drho + 0.5);

  // ------------------------------------------------------------------
  // setup frho arrays
  // ------------------------------------------------------------------
  frho.resize(eam_paras.nrho + 1);

  rbmd::Real r, p, cof1, cof2, cof3, cof4;
  for (m = 1; m <= eam_paras.nrho; m++)
  {
    r = (m - 1) * eam_paras.drho;
    p = r / file.drho + 1.0;
    k = static_cast<int>(p);
    k = MIN(k, file.nrho - 2);
    k = MAX(k, 2);
    p -= k;
    p = MIN(p, 2.0);
    cof1 = -sixth * p * (p - 1.0) * (p - 2.0);
    cof2 = 0.5 * (p * p - 1.0) * (p - 2.0);
    cof3 = -0.5 * p * (p + 1.0) * (p - 2.0);
    cof4 = sixth * p * (p * p - 1.0);
    frho[m] = cof1 * file.frho[k - 1] + cof2 * file.frho[k] + cof3 * file.frho[k + 1] +
      cof4 * file.frho[k + 2];
  }

  // ------------------------------------------------------------------
  // setup rhor arrays
  // ------------------------------------------------------------------
  rhor.resize(eam_paras.nrho + 1);
  for (m = 1; m <= eam_paras.nr; m++)
  {
    r = (m - 1) * eam_paras.dr;
    p = r / file.dr + 1.0;
    k = static_cast<int>(p);
    k = MIN(k, file.nr - 2);
    k = MAX(k, 2);
    p -= k;
    p = MIN(p, 2.0);
    auto cof1 = -sixth * p * (p - 1.0) * (p - 2.0);
    auto cof2 = 0.5 * (p * p - 1.0) * (p - 2.0);
    auto cof3 = -0.5 * p * (p + 1.0) * (p - 2.0);
    auto cof4 = sixth * p * (p * p - 1.0);
    rhor[m] = cof1 * file.rhor[k - 1] + cof2 * file.rhor[k] + cof3 * file.rhor[k + 1] +
      cof4 * file.rhor[k + 2];
  }

  // ------------------------------------------------------------------
  // setup z2r arrays
  // ------------------------------------------------------------------
  z2r.resize(eam_paras.nr + 1);

  double zri;
  for (m = 1; m <= eam_paras.nr; m++)
  {
    r = (m - 1) * eam_paras.dr;

    p = r / file.dr + 1.0;
    k = static_cast<int>(p);
    k = MIN(k, file.nr - 2);
    k = MAX(k, 2);
    p -= k;
    p = MIN(p, 2.0);
    cof1 = -sixth * p * (p - 1.0) * (p - 2.0);
    cof2 = 0.5 * (p * p - 1.0) * (p - 2.0);
    cof3 = -0.5 * p * (p + 1.0) * (p - 2.0);
    cof4 = sixth * p * (p * p - 1.0);
    zri = cof1 * file.zr[k - 1] + cof2 * file.zr[k] + cof3 * file.zr[k + 1] + cof4 * file.zr[k + 2];

    z2r[m] = 27.2 * 0.529 * zri * zri;
  }
}

void EAM::array2spline() {
  _h_frho_spline.resize(eam_paras.nrho + 1);
  _h_rhor_spline.resize(eam_paras.nrho + 1);
  _h_z2r_spline.resize(eam_paras.nr + 1);

  interpolate(eam_paras.nrho, eam_paras.drho, frho, _h_frho_spline);
  interpolate(eam_paras.nr, eam_paras.dr, rhor, _h_rhor_spline);
  interpolate(eam_paras.nr, eam_paras.dr, z2r, _h_z2r_spline);
}

void EAM::interpolate(rbmd::Id n, rbmd::Real delta, std::vector<rbmd::Real>& f,
                       thrust::host_vector<Real7>& spline) {

  for (int m = 1; m <= n; m++)
  {
    spline[m][6] = f[m];  //f(x)
  }

  spline[1][5] = spline[2][6] -
    spline[1][6]; //f'(x) = (f(x + h) - f(x)) / h    [5] is the coefficient of the first derivative(energy expression)
  spline[2][5] = 0.5 * (spline[3][6] - spline[1][6]);
  spline[n - 1][5] = 0.5 * (spline[n][6] - spline[n - 2][6]);
  spline[n][5] = spline[n][6] - spline[n - 1][6];

  for (int m = 3; m <= n - 2; m++)
  {
    spline[m][5] =
      ((spline[m - 2][6] - spline[m + 2][6]) + 8.0 * (spline[m + 1][6] - spline[m - 1][6])) /
      12.0; //further sample points for a more accurate estimate
  }

  for (int m = 1; m <= n - 1; m++)
  {
    spline[m][4] = 3.0 * (spline[m + 1][6] - spline[m][6]) - 2.0 * spline[m][5] -
      spline[m + 1][5]; //[4] is the coefficient of the second derivative
    spline[m][3] = spline[m][5] + spline[m + 1][5] -
      2.0 * (spline[m + 1][6] - spline[m][6]); // [3] is the coefficient of the third derivative
  }

  spline[n][4] = 0.0;
  spline[n][3] = 0.0; //The second and third derivative coefficients at the last sample point are zero,
  // To make the interpolation curve smoother at both ends, the higher derivative coefficient at the boundary can be set to zero.
  // This is because spline interpolation typically uses higher-order polynomial interpolation
  // at the inner sample points and lower-order polynomials at the boundaries to ensure smoothness.

  for (int m = 1; m <= n; m++)
  {
    spline[m][2] = spline[m][5] / delta;       //The coefficient of the second derivative(force expression)
    spline[m][1] = 2.0 * spline[m][4] / delta; //The coefficient of the first derivative
    spline[m][0] = 3.0 * spline[m][3] / delta; //The coefficient of the zero derivative (i.e. the value of the function).
  }
}

void EAM::SetEAM() {

  _d_frho_spline.resize(eam_paras.nrho + 1);
  _d_rhor_spline.resize(eam_paras.nrho + 1);
  _d_z2r_spline.resize(eam_paras.nr + 1);

  //H2D
  thrust::copy(_h_frho_spline.begin(),
  _h_frho_spline.end(), _d_frho_spline.begin());

  thrust::copy(_h_rhor_spline.begin(),
_h_rhor_spline.end(), _d_rhor_spline.begin());

  thrust::copy(_h_z2r_spline.begin(),
_h_z2r_spline.end(), _d_z2r_spline.begin());
}

void EAM::EAMVerlet() {
    //neighbor_list_build
  auto start = std::chrono::high_resolution_clock::now();
  _list = _neighbor_list_builder->Build();

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration = end - start;
  TimingStatistics::Instance().record("Neighbor-List",duration.count());

  //EAM_fp
  auto num_atoms = *(_structure_info_data->_num_atoms);
  thrust::device_vector<rbmd::Real> eam_rho(num_atoms);
  thrust::device_vector<rbmd::Real> eam_fp(num_atoms);

  thrust::device_vector<rbmd::Real> d_energy_embedding(1, 0.0);
  thrust::device_vector<rbmd::Real> d_energy_pair(1, 0.0);

  auto start_f = std::chrono::high_resolution_clock::now();
  op::ComputeEAMForceVerlet<device::DEVICE_GPU>()(
*_box, eam_paras ,file.cut_off, num_atoms,
thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
thrust::raw_pointer_cast(_device_data->_d_atoms_id.data()),
thrust::raw_pointer_cast(_list->_start_idx.data()),
thrust::raw_pointer_cast(_list->_end_idx.data()),
thrust::raw_pointer_cast(_list->_d_neighbors.data()),
thrust::raw_pointer_cast(_d_rhor_spline.data()),
thrust::raw_pointer_cast(_d_frho_spline.data()),
thrust::raw_pointer_cast(_d_z2r_spline.data()),
thrust::raw_pointer_cast(_device_data->_d_px.data()),
thrust::raw_pointer_cast(_device_data->_d_py.data()),
thrust::raw_pointer_cast(_device_data->_d_pz.data()),
thrust::raw_pointer_cast(eam_rho.data()),
thrust::raw_pointer_cast(eam_fp.data()),
thrust::raw_pointer_cast(_device_data->_d_fx.data()),
thrust::raw_pointer_cast(_device_data->_d_fy.data()),
thrust::raw_pointer_cast(_device_data->_d_fz.data()),
thrust::raw_pointer_cast(d_energy_embedding.data()),
thrust::raw_pointer_cast(d_energy_pair.data()));

  // D2H
  thrust::host_vector<rbmd::Real> h_energy_embedding(d_energy_embedding);
  thrust::host_vector<rbmd::Real> h_energy_pair(d_energy_pair);
  _e_embedding = h_energy_embedding[0] / num_atoms;
  _e_pair = ( h_energy_pair[0]) / num_atoms;

  auto end_f = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration_f = end_f - start_f;
  TimingStatistics::Instance().record("Short-Range",duration_f.count());
}

void EAM::EAMRBL() {
  //neighbor_list_build
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
  thrust::device_vector<rbmd::Real> eam_rho(num_atoms);
  thrust::device_vector<rbmd::Real> eam_fp(num_atoms);

 //fp
  op::ComputeFpRBL<device::DEVICE_GPU>()(
*_box, eam_paras ,r_core ,file.cut_off, num_atoms,neighbor_sample_num,
_rbl_list->_selection_frequency,
thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
thrust::raw_pointer_cast(_device_data->_d_atoms_id.data()),
thrust::raw_pointer_cast(_rbl_list->_start_idx.data()),
   thrust::raw_pointer_cast(_rbl_list->_end_idx.data()),
   thrust::raw_pointer_cast(_rbl_list->_d_neighbors.data()),
   thrust::raw_pointer_cast(_rbl_list->_d_random_neighbor.data()),
   thrust::raw_pointer_cast(_rbl_list->_d_random_neighbor_num.data()),
thrust::raw_pointer_cast(_d_rhor_spline.data()),
thrust::raw_pointer_cast(_d_frho_spline.data()),
thrust::raw_pointer_cast(_device_data->_d_px.data()),
thrust::raw_pointer_cast(_device_data->_d_py.data()),
thrust::raw_pointer_cast(_device_data->_d_pz.data()),
thrust::raw_pointer_cast(eam_fp.data()));

  //EAMForce
  op::ComputeEAMForceRBL<device::DEVICE_GPU>()(
*_box, eam_paras ,r_core,file.cut_off, num_atoms,neighbor_sample_num,
_rbl_list->_selection_frequency,
thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
thrust::raw_pointer_cast(_device_data->_d_atoms_id.data()),
thrust::raw_pointer_cast(_rbl_list->_start_idx.data()),
thrust::raw_pointer_cast(_rbl_list->_end_idx.data()),
thrust::raw_pointer_cast(_rbl_list->_d_neighbors.data()),
thrust::raw_pointer_cast(_rbl_list->_d_random_neighbor.data()),
thrust::raw_pointer_cast(_rbl_list->_d_random_neighbor_num.data()),
thrust::raw_pointer_cast(_d_rhor_spline.data()),
thrust::raw_pointer_cast(_d_z2r_spline.data()),
thrust::raw_pointer_cast(_device_data->_d_px.data()),
thrust::raw_pointer_cast(_device_data->_d_py.data()),
thrust::raw_pointer_cast(_device_data->_d_pz.data()),
thrust::raw_pointer_cast(eam_fp.data()),
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
  if ("yes" == _energy_rbl_flag ) {
    ComputEAMEnergy();
  }

}

void EAM::ComputEAMEnergy() {
  //neighbor_list_build
  _list = _neighbor_list_builder->Build();


  //
  auto num_atoms = *(_structure_info_data->_num_atoms);
  thrust::device_vector<rbmd::Real> eam_rho(num_atoms);
  thrust::device_vector<rbmd::Real> eam_fp(num_atoms);
  thrust::device_vector<rbmd::Real> d_energy_embedding(1, 0.0);
  thrust::device_vector<rbmd::Real> d_energy_pair(1, 0.0);

  op::ComputeEAMEnergy<device::DEVICE_GPU>()(
*_box, eam_paras ,file.cut_off, num_atoms,
thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
thrust::raw_pointer_cast(_device_data->_d_atoms_id.data()),
thrust::raw_pointer_cast(_list->_start_idx.data()),
thrust::raw_pointer_cast(_list->_end_idx.data()),
thrust::raw_pointer_cast(_list->_d_neighbors.data()),
thrust::raw_pointer_cast(_d_rhor_spline.data()),
thrust::raw_pointer_cast(_d_frho_spline.data()),
thrust::raw_pointer_cast(_d_z2r_spline.data()),
thrust::raw_pointer_cast(_device_data->_d_px.data()),
thrust::raw_pointer_cast(_device_data->_d_py.data()),
thrust::raw_pointer_cast(_device_data->_d_pz.data()),
thrust::raw_pointer_cast(eam_rho.data()),
thrust::raw_pointer_cast(eam_fp.data()),
thrust::raw_pointer_cast(d_energy_embedding.data()),
thrust::raw_pointer_cast(d_energy_pair.data()));

  // D2H
  thrust::host_vector<rbmd::Real> h_energy_embedding(d_energy_embedding);
  thrust::host_vector<rbmd::Real> h_energy_pair(d_energy_pair);
  _e_embedding = h_energy_embedding[0] / num_atoms;
  _e_pair = ( h_energy_pair[0]) / num_atoms;

}

void EAM::SumForces() {
  if ("RBL" ==_neighbor_type) {
    EAMRBL();
  }
  else {
    EAMVerlet();
  }

  //
//   auto atom_id_to_idx =
// LinkedCellLocator::GetInstance().GetLinkedCell()->_atom_id_to_idx;
//    thrust::host_vector<rbmd::Real> h_fx(num_atoms);
//    thrust::host_vector<rbmd::Real> h_fy(num_atoms);
//    thrust::host_vector<rbmd::Real> h_fz(num_atoms);
//     h_fx = _device_data->_d_fx;
//     h_fy = _device_data->_d_fy;
//     h_fz = _device_data->_d_fz;
//
//    std::ofstream fx("fx.txt");
//    if (fx.is_open()) {
//      for (rbmd::Id i = 0; i < num_atoms; ++i) {
//        auto idx = atom_id_to_idx[i];
//        fx << i  << " " << h_fx[idx] << " " << h_fy[idx]
//          << " " << h_fz[idx]  << "\n";
//      }
//      fx.close();
//   }


}

void EAM::EvaluatePotentialenergy() {

  _e_pe = _e_embedding + _e_pair;
  ThermoStats::Instance().AddThermoData("total-potential-energy",_e_pe);

  //out
  auto interval = DataManager::getInstance().getConfigData()->Get<rbmd::Id>(
"interval", "outputs", "thermo_out");

  std::ofstream outfile("thermo.txt", std::ios::app);
  if (outfile.tellp() == 0) {
    outfile << "step  e_pe" << std::endl;
  }
  if (test_current_step % interval == 0) {
    outfile << test_current_step << " " <<  _e_pe << std::endl;
  }
  outfile.close();
}


