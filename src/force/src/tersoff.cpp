#include "tersoff.h"

#include <regex>
#include <unordered_map>
#include <thrust/device_ptr.h>
#include "thrust/sort.h"

#include "../../common/device_types.h"
#include "../../common/rbmd_define.h"
#include "../../common/types.h"
#include "../common/unit_factor.h"
#include "tersoff_op/tersoff_op.h"

#include "neighbor_list/include/linked_cell/linked_cell_locator.h"
#include "neighbor_list/include/neighbor_list_builder/half_neighbor_list_builder.h"
#include "neighbor_list/include/neighbor_list_builder/full_neighbor_list_builder.h"
#include "neighbor_list/include/neighbor_list_builder/rbl_full_neighbor_list_builder.h"

#include "common/thermo_stats.hpp"
#include "common/timing_statistics.hpp"

extern int test_current_step;
extern std::map<std::string, UNIT> unit_factor_map;

TerSoff::TerSoff():
  _elements(nullptr)
{
  _rbl_neighbor_list_builder = std::make_shared<RblFullNeighborListBuilder>();
  _neighbor_list_builder = std::make_shared<FullNeighborListBuilder>();
  _h_params = nullptr;

  const auto& config = DataManager::getInstance().getConfigData();
  if (config->PathExists({"hyper_parameters", "neighbor","interval" }))
  {
    _interval = config->Get<rbmd::Id>(
  "interval", "hyper_parameters", "neighbor");
  }

  //map resize
  auto atoms_type = *(_structure_info_data->_num_atoms_type);
  _map.resize(atoms_type +1 ,-1);
  _setflag.resize(atoms_type + 1,std::vector<rbmd::Id>(atoms_type + 1, 0));
  std::remove("thermo.txt");
}

TerSoff::~TerSoff()
{
  free(_h_params);
  FREE(d_params);

  if (_elements) {
    for (int i = 0; i < _nelements; i++) {delete[] _elements[i];}
  }
  delete[] _elements;

}

void TerSoff::Init()
{
  const auto& config = DataManager::getInstance().getConfigData();
  auto shift_flag =config->PathExists({"hyper_parameters", "shift_value" });
  if (shift_flag) {
    _shift.shift_flag =1;
    _shift.shift_value= config->Get<rbmd::Real>("shift_value", "hyper_parameters");
  }

  _cut_off = config->Get<rbmd::Real>("cut_off", "hyper_parameters", "neighbor");
  //
  int narg = 0;
  char** arg = nullptr;
  bool update_setflag = true;
  auto potential_elements = DataManager::getInstance().getConfigData()->Get
    <std::string>("potential_elements", "hyper_parameters", "force_field");
  ReadPotentialElements(potential_elements, narg, &arg);
  Element2Type(narg-2,arg+2,update_setflag);
  //
  std::cout << "narg: " << narg << std::endl;
  for (int i = 0; i < narg; ++i) {
    std::cout << "arg[" << i << "]: " << arg[i] << std::endl;
  }


  //
  auto potential_file = DataManager::getInstance().getConfigData()->Get
    <std::string>("potential_file", "hyper_parameters", "force_field");
  _potential_file.open(potential_file);

  if (!_potential_file.is_open()) {
    Logger::Instance().error("\033[31m Failed to open potential file: {}\033[0m", potential_file );
    exit(EXIT_FAILURE); //
  } else {
    Logger::Instance().info(" Successfully opened potential file:  {}", potential_file);
  }

  ReadPotentialFile_fix(_potential_file);
  _potential_file.close();

  SetupParams();

  //
  CHECK_RUNTIME(MALLOC(&d_params, _nparams * sizeof(TersoffParams)));
  MEMCPY(d_params,_h_params,_nparams * sizeof(TersoffParams),H2D);


  //
  auto start_list = std::chrono::high_resolution_clock::now();
  if (test_current_step == 0) {
    _list = _neighbor_list_builder->Build(_cutmax);
  }
  auto end_list = std::chrono::high_resolution_clock::now();
  _duration_list_init = end_list - start_list;

}

void TerSoff::Execute()
{
  ComputeTersoff();

  EvaluatePotentialenergy();
}

void TerSoff::ComputeTersoff() {
  //neighbor_list_build
  auto start = std::chrono::high_resolution_clock::now();
  if (test_current_step>0) {
    if (test_current_step  % _interval == 0) {
      _list = _neighbor_list_builder->Build(_cutmax);
    }
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration = end - start;
  duration = _duration_list_init + duration;
  TimingStatistics::Instance().record("Neighbor-List",duration.count());

  //Tersoff
  auto start_f = std::chrono::high_resolution_clock::now();

  auto num_atoms = *(_structure_info_data->_num_atoms);
  thrust::device_vector<rbmd::Real> d_total_energy(1, 0.0);
  auto atom_id_to_idx =
    LinkedCellLocator::GetInstance().GetLinkedCell()->_atom_id_to_idx;

  thrust::fill(_device_data->_d_fx.begin(),
  _device_data->_d_fx.end(), 0.0f);
  thrust::fill(_device_data->_d_fy.begin(),
    _device_data->_d_fy.end(), 0.0f);
  thrust::fill(_device_data->_d_fz.begin(),
    _device_data->_d_fz.end(), 0.0f);

  op::TerSoff<device::DEVICE_GPU>()(
    *_box,d_params,_shift,_cutmax,num_atoms,_nelements,
    thrust::raw_pointer_cast(atom_id_to_idx.data()),
    thrust::raw_pointer_cast(_device_data->_d_atoms_id.data()),
    thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
    thrust::raw_pointer_cast(_d_map.data()),
    thrust::raw_pointer_cast(_d_elem3param.data()),
    thrust::raw_pointer_cast(_list->_start_idx.data()),
    thrust::raw_pointer_cast(_list->_end_idx.data()),
    thrust::raw_pointer_cast(_list->_d_neighbors.data()),
    thrust::raw_pointer_cast(_device_data->_d_px.data()),
thrust::raw_pointer_cast(_device_data->_d_py.data()),
thrust::raw_pointer_cast(_device_data->_d_pz.data()),
thrust::raw_pointer_cast(_device_data->_d_fx.data()),
thrust::raw_pointer_cast(_device_data->_d_fy.data()),
thrust::raw_pointer_cast(_device_data->_d_fz.data()),
thrust::raw_pointer_cast(_device_data->_d_flat_virial_lj.data()),
thrust::raw_pointer_cast(d_total_energy.data()));

  auto end_f = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration_f = end_f - start_f;
  TimingStatistics::Instance().record("Short-Range",duration_f.count());
  // D2H
  thrust::host_vector<rbmd::Real> h_total_evdwl(d_total_energy);
  _e_vdwl = h_total_evdwl[0]/num_atoms;

  //sum virial_special_lj on host
  ReduceVirial(num_atoms,_device_data->_d_flat_virial_lj,
_device_data->_d_virial_lj);

  // thrust::host_vector<rbmd::Real> h_fx;
  // thrust::host_vector<rbmd::Real> h_fy;
  // thrust::host_vector<rbmd::Real> h_fz;
  // h_fx.reserve(num_atoms);
  // h_fy.reserve(num_atoms);
  // h_fz.reserve(num_atoms);
  // h_fx = _device_data->_d_fx;
  // h_fy = _device_data->_d_fy;
  // h_fz = _device_data->_d_fz;
  //
  // std::ofstream force("force.txt");
  // if (force.is_open()) {
  //   for (rbmd::Id i = 0; i <h_fx.size(); ++i) {
  //     auto idx = atom_id_to_idx[i];
  //     force << i  << " " << h_fx[idx]  << " "<< h_fy[idx]  << " " << h_fz[idx] << "\n";
  //   }
  //   force.close();
  // }

}

void TerSoff::ReadPotentialElements(const std::string& potential_elements,
  int& narg, char*** arg){

  //
  std::istringstream iss(potential_elements);
  std::vector<std::string> tokens;
  std::string token;
  while (iss >> token) {
    tokens.push_back(token);
  }

  //
  narg = tokens.size();
  if (narg == 0) {
    *arg = nullptr;
    return;
  }

  //
  *arg = new char*[narg];
  for (int i = 0; i < narg; ++i) {
    (*arg)[i] = new char[tokens[i].size() + 1];
    std::strcpy((*arg)[i], tokens[i].c_str());
  }
}

void TerSoff::Element2Type(rbmd::Id narg, char **arg, bool update_setflag)
{
  rbmd::Id i, j;
  const rbmd::Id ntypes = *(_structure_info_data->_num_atoms_type);

  if (narg != ntypes)
  {
    Logger::Instance().error("\033[31m Incorrect element mapping for "
        "tersoff coefficients. The total number of atom types is {}\033[0m", ntypes );
    exit(EXIT_FAILURE); //
  }

  // //
  if (_elements) {
    for (i = 0; i < _nelements; i++) {
      delete[] _elements[i];
    }
    delete[] _elements;
  }

  //
  _elements = new char*[ntypes];
  for (i = 0; i < ntypes; i++) {
    _elements[i] = nullptr;
  }

  _nelements = 0;
  _map[0] = -1;  //

  for (i = 1; i <= narg; i++) {
    std::string entry = arg[i-1];
    if (entry == "NULL") {
      _map[i] = -1;
      continue;
    }
    for (j = 0; j < _nelements; j++) {
      //if (strcmp(entry.c_str(), _elements[j]) == 0)
      if (entry == _elements[j])
        break;
    }
    _map[i] = j;

    if (j == _nelements) {
      _elements[j] = StrDup(entry);
      _nelements++;
    }
  }
  // for (int i = 0; i < _map.size(); ++i) {
  //     std::cout << "map: "   <<  i   << "  "<< _map[i] << " " << std::endl;
  // }
  // for (i = 0; i < ntypes; i++) {
  //  std::cout <<"elements: " << i << " "<< _elements[i]  << std::endl;
  // }

  _d_map = _map;
  //
  update_setflag = true;
  if (update_setflag) {
    int count = 0;
    for (i = 1; i <= ntypes; i++) {
      for (j = i; j <= ntypes; j++) {
        _setflag[i][j] = 0;
        if ((_map[i] >= 0) && (_map[j] >= 0)) {
          _setflag[i][j] = 1;
          count++;
        }
      }
    }
    if (count == 0) {
      Logger::Instance().error("\033[31m Incorrect args for tersoff coefficients\033[0m" );
    }
  }
}

void TerSoff::ReadPotentialFile_fix(std::ifstream& file)
{
    _nparams = _maxparam = 0;
    std::string accumulatedLine;  //

    std::string line;
    while (std::getline(file, line))
    {
        // trim
        auto start = line.find_first_not_of(" \t");
        if (start == std::string::npos) {
            continue;
        }
        auto end = line.find_last_not_of(" \t");
        line = line.substr(start, end - start + 1);

        //
        if (line[0] == '#') {
            continue;
        }
        //
        if (accumulatedLine.empty()) {
            accumulatedLine = line;
        } else {
            accumulatedLine += " " + line;
        }
        //
        bool parseSuccess = false;
        try {
            std::istringstream iss(accumulatedLine);

            //
            std::string iname, jname, kname;
          if (!(iss >> iname >> jname >> kname))
             continue;
            //
            int ielement = -1, jelement = -1, kelement = -1;
            for (int i = 0; i < _nelements; ++i) {
                if (iname == _elements[i]) ielement = i;
                if (jname == _elements[i]) jelement = i;
                if (kname == _elements[i]) kelement = i;
            }

          //
          if (ielement == -1 || jelement == -1 || kelement == -1) {
              //
              accumulatedLine.clear();
              continue;
          }

            //
            rbmd::Id DELTA = 4;
            if (_nparams == _maxparam) {
                _maxparam += DELTA;
                _h_params = (TersoffParams *)realloc(_h_params, _maxparam * sizeof(TersoffParams));
                memset(_h_params + _nparams, 0, DELTA * sizeof(TersoffParams));
            }

            //
            _h_params[_nparams].ielement = ielement;
            _h_params[_nparams].jelement = jelement;
            _h_params[_nparams].kelement = kelement;

            //
            if (!(iss >> _h_params[_nparams].m >> _h_params[_nparams].gamma >> _h_params[_nparams].lambda3
                  >> _h_params[_nparams].c >> _h_params[_nparams].d >> _h_params[_nparams].costheta0
                  >> _h_params[_nparams].n >> _h_params[_nparams].beta >> _h_params[_nparams].lambda2
                  >> _h_params[_nparams].B >> _h_params[_nparams].R >> _h_params[_nparams].D
                  >> _h_params[_nparams].lambda1 >> _h_params[_nparams].A))
              continue;
            //
            _h_params[_nparams].m_int = rbmd::Id(_h_params[_nparams].m);

            // std::cout << "Read parameters: "
            //           << "m = " << _h_params[_nparams].m << ", "
            //           << "gamma = " << _h_params[_nparams].gamma << ", "
            //           << "lambda3 = " << _h_params[_nparams].lambda3 << ", "
            //           << "c = " << _h_params[_nparams].c << ", "
            //           << "d = " << _h_params[_nparams].d << ", "
            //           << "costheta0 = " << _h_params[_nparams].costheta0 << ", "
            //           << "n = " << _h_params[_nparams].n << ", "
            //           << "beta = " << _h_params[_nparams].beta << ", "
            //           << "lambda2 = " << _h_params[_nparams].lambda2 << ", "
            //           << "B = " << _h_params[_nparams].B << ", "
            //           << "R = " << _h_params[_nparams].R << ", "
            //           << "D = " << _h_params[_nparams].D << ", "
            //           << "lambda1 = " << _h_params[_nparams].lambda1 << ", "
            //           << "A = " << _h_params[_nparams].A
            //           << std::endl;

            //
            bool unit_convert_flag = true;
            rbmd::Real conversion_factor = 1.0;
            if (unit_convert_flag) {
                _h_params[_nparams].A *= conversion_factor;
                _h_params[_nparams].B *= conversion_factor;
            }

            //
            parseSuccess = true;
            accumulatedLine.clear();
           //++_nparams;

        } catch (const std::exception& e) {
            std::cerr << "Error parsing accumulated line: " << accumulatedLine << "\n"
                      << "Error: " << e.what() << std::endl;
            accumulatedLine.clear(); //
        }
        _nparams++;
        //
        if (!parseSuccess && accumulatedLine.length() > 1000) {
            accumulatedLine.clear();
        }
    }
    //
    if (_nparams == 0) {
        Logger::Instance().error("\033[31m No valid parameters found in the potential file\033[0m");
        exit(EXIT_FAILURE);
    }
    std::cout << "_nparams: " << _nparams << std::endl;
}

void TerSoff::ReadPotentialFile(std::ifstream& file)
{

    _nparams = _maxparam = 0;

    std::string line;
    while (std::getline(file, line))
    {
      //
      line = line.substr(line.find_first_not_of(" \t"), line.find_last_not_of(" \t") + 1);

      //
      if (line.empty() || line[0] == '#'){
        continue;
      }

      std::istringstream iss(line);
      try {
            //
            std::string iname, jname, kname;
            iss >> iname >> jname >> kname;

            //
            int ielement, jelement, kelement;
            for (ielement = 0; ielement < _nelements; ++ielement)
                if (iname == _elements[ielement]) break;
            if (ielement == _nelements) continue;

            for (jelement = 0; jelement < _nelements; ++jelement)
                if (jname == _elements[jelement]) break;
            if (jelement == _nelements) continue;

            for (kelement = 0; kelement < _nelements; ++kelement)
                if (kname == _elements[kelement]) break;
            if (kelement == _nelements) continue;

            rbmd::Id  DELTA =4;
            if (_nparams == _maxparam)
            {
              _maxparam += DELTA;
              _h_params = (TersoffParams *) realloc(_h_params,_maxparam*sizeof(TersoffParams));
              memset(_h_params + _nparams, 0, DELTA*sizeof(TersoffParams));
           }

            //
            _h_params[_nparams].ielement = ielement;
            _h_params[_nparams].jelement = jelement;
            _h_params[_nparams].kelement = kelement;

            //
            iss >> _h_params[_nparams].m >> _h_params[_nparams].gamma >> _h_params[_nparams].lambda3
                >> _h_params[_nparams].c >> _h_params[_nparams].d>> _h_params[_nparams].costheta0 >>
          _h_params[_nparams].n >>_h_params[_nparams].beta>> _h_params[_nparams].lambda2>>
          _h_params[_nparams].B >> _h_params[_nparams].R >>_h_params[_nparams].D >>
          _h_params[_nparams].lambda1 >> _h_params[_nparams].A;

          _h_params[_nparams].m_int = rbmd::Id(_h_params[_nparams].m);//

            std::cout << "Read parameters: "
            << "m = " << _h_params[_nparams].m << ", "
            << "gamma = " << _h_params[_nparams].gamma << ", "
            << "lambda3 = " << _h_params[_nparams].lambda3 << ", "
            << "c = " << _h_params[_nparams].c << ", "
            << "d = " << _h_params[_nparams].d << ", "
            << "costheta0 = " << _h_params[_nparams].costheta0 << ", "
            << "n = " << _h_params[_nparams].n << ", "
            << "beta = " << _h_params[_nparams].beta << ", "
            << "lambda2 = " << _h_params[_nparams].lambda2 << ", "
            << "B = " << _h_params[_nparams].B << ", "
            << "R = " << _h_params[_nparams].R << ", "
            << "D = " << _h_params[_nparams].D << ", "
            << "lambda1 = " << _h_params[_nparams].lambda1 << ", "
            << "A = " << _h_params[_nparams].A
            << std::endl;

            //
            bool unit_convert_flag =true;
            rbmd::Real conversion_factor=1.0;
            if (unit_convert_flag) {
                _h_params[_nparams].A *= conversion_factor;
                _h_params[_nparams].B *= conversion_factor;
            }

            //
         } catch (const std::exception& e) {
           std::cerr << "Error parsing line [" << _nparams << "]: " << line << "\n"
                        << "Error: " << e.what() << std::endl;
         }
         ++_nparams;
    }
    //
    if (_nparams == 0) {
      Logger::Instance().error("\033[31m No valid parameters found in the potential file\033[0m" );
      exit(EXIT_FAILURE); //
    }

  //std::cout<< "_nparams: " <<_nparams<<std::endl;
}

void TerSoff::SetupParams()
{
    int i, j, k, m, n;
    //
  _elem3param = std::vector<std::vector<std::vector<rbmd::Id>>>(
      _nelements,
      std::vector<std::vector<rbmd::Id>>(
          _nelements,
          std::vector<rbmd::Id>(_nelements, -1)
      ));

    //
    for (i = 0; i < _nelements; i++) {
        for (j = 0; j < _nelements; j++) {
            for (k = 0; k < _nelements; k++) {
                n = -1; //
                for (m = 0; m < _nparams; m++) {
                    if (_h_params[m].ielement == i &&
                        _h_params[m].jelement == j &&
                        _h_params[m].kelement == k) {
                        if (n >= 0) {
                        throw std::runtime_error(
                            std::string("Duplicate entry in potential file for _elements: ") +
                            _elements[i] + " " + _elements[j] + " " + _elements[k]);
                        }
                        n = m; //
                    }
                }
                if (n < 0) {
                  throw std::runtime_error(
                    std::string("Missing entry in potential file for _elements: ") +
                    _elements[i] + " " + _elements[j] + " " + _elements[k]);
                }
                _elem3param[i][j][k] = n; //
            }
        }
    }

  // for (int i = 0; i <_nelements; i++) {
  //   for (int j = 0; j < _nelements; j++) {
  //     for (int k = 0; k < _nelements; k++) {
  //       std::cout << "_elem3param[" << i << "][" << j << "][" << k << "] = "
  //                 << _elem3param[i][j][k] << std::endl;
  //     }
  //   }
  // }

  //h_elem3param_flat
  std::vector<rbmd::Id> _h_elem3param_flat(_nelements * _nelements * _nelements, -1);
  for (int i = 0; i < _nelements; i++) {
    for (int j = 0; j < _nelements; j++) {
      for (int k = 0; k < _nelements; k++) {
        int flat_idx = i * _nelements * _nelements + j * _nelements + k;
        _h_elem3param_flat[flat_idx] = _elem3param[i][j][k];
      }
    }
  }

   _d_elem3param = thrust::device_vector<rbmd::Id>(_h_elem3param_flat.begin(), _h_elem3param_flat.end());


  //c1  c2 c3  c4
  for (int i = 0; i < _nparams; i++) {
    _h_params[i].cut = _h_params[i].R + _h_params[i].D;
    _h_params[i].cutsq = _h_params[i].cut * _h_params[i].cut;

    if (_h_params[i].n > 0.0) {
      _h_params[i].c1 = POW(2.0 * _h_params[i].n * 1.0e-16, -1.0 / _h_params[i].n);
      _h_params[i].c2 = POW(2.0 * _h_params[i].n * 1.0e-8, -1.0 / _h_params[i].n);
      _h_params[i].c3 = 1.0 / _h_params[i].c2;
      _h_params[i].c4 = 1.0 / _h_params[i].c1;
    } else {
      _h_params[i].c1 = _h_params[i].c2 = _h_params[i].c3 = _h_params[i].c4 = 0.0;
    }
  }

  //cutmax
  _cutmax = 0.0;
  for (int i = 0; i < _nparams; i++) {
    if (_h_params[i].cut > _cutmax) _cutmax = _h_params[i].cut;
  }
  Logger::Instance().info(" max cut_off= {}", _cutmax);
 }

void TerSoff::EvaluatePotentialenergy()
{
  _e_pe = _e_vdwl;

  ThermoStats::Instance().AddThermoData("total-potential-energy",_e_pe);

  //out
  auto interval = DataManager::getInstance().getConfigData()->Get<rbmd::Id>(
"interval", "outputs", "thermo_out");

  std::ofstream outfile("thermo.txt", std::ios::app);
  if (outfile.tellp() == 0) {
    outfile << "step  e_vdwl  e_pe" << std::endl;
  }
  if (test_current_step % interval == 0) {
    outfile << test_current_step << " " << _e_vdwl  << " " <<  _e_pe << std::endl;
  }
  outfile.close();
}