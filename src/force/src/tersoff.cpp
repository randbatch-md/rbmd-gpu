#include "tersoff.h"
#include <sstream> 
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
// #include <hipcub/hipcub.hpp>
// #include <hipcub/backend/rocprim/block/block_reduce.hpp>


extern int test_current_step;
extern std::map<std::string, UNIT> unit_factor_map;

TerSoff::TerSoff():
  _elements(nullptr)
{
  _rbl_neighbor_list_builder = std::make_shared<RblFullNeighborListBuilder>();
  _neighbor_list_builder = std::make_shared<FullNeighborListBuilder>();
  _params = nullptr;

  //
  _map.resize(100,-1);

  auto atoms_type = *(_structure_info_data->_num_atoms_type);
  _setflag.resize(atoms_type + 1,std::vector<rbmd::Id>(atoms_type + 1, 0));

}

TerSoff::~TerSoff()
{
  free(_params);

  if (_elements) {
    for (int i = 0; i < _nelements; i++) {delete[] _elements[i];}
  }
  delete[] _elements;

}

void TerSoff::Init()
{
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
    std::cerr << "Error: Failed to open potential file: " << potential_file << std::endl;
    return;
  } else {
    std::cout << "Successfully opened potential file: " << potential_file << std::endl;
  }

  ReadPotentialFile_fix(_potential_file);
  _potential_file.close();

  SetupParams();
}

void TerSoff::Execute()
{
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
     std::cerr <<
       "Incorrect element mapping for tersoff coefficients"
    << std::endl;;
  }

  // unordered_map
  std::unordered_map<std::string, int> element_to_index;

  //
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
    for (j = 0; j < _nelements; j++)
      if (entry == _elements[j]) break;
    _map[i] = j;
    if (j == _nelements) {
      _elements[j] = StrDup(entry);
      _nelements++;
    }
  }


  // //
  // for (i = 1; i <= narg; i++)
  // {
  //   std::string entry = arg[i - 1];
  //   if (entry == "NULL") {
  //     _map[i] = -1;  // "NULL" 映射为 -1
  //     continue;
  //   }
  //
  //   //
  //   auto it = element_to_index.find(entry);
  //   if (it != element_to_index.end()) {
  //     _map[i] = it->second;
  //   } else {
  //     //
  //     _elements[_nelements] = StrDup(entry);
  //     element_to_index[entry] = _nelements;
  //     _map[i] = _nelements;
  //     _nelements++;
  //   }
  // }

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
      std::cerr << "Incorrect args for tersoff coefficients"<< std::endl;
    }
  }
}

// void TerSoff::ReadPotentialFile(std::ifstream& file)
// {
//   if (!file.is_open())
//   {
//     std::cerr << "Unable to open the file." << std::endl;
//   }
//
//   TersoffData tersoffData;
//   std::string line;
//   std::regex dataRegex(R"((\S+)\s+(\S+)\s+(\S+)\s+([\S\s]+))");
//
//   //
//   while (std::getline(file, line))
//   {
//     //
//     line = line.substr(line.find_first_not_of(" \t"), line.find_last_not_of(" \t") + 1);
//
//     //
//     if (line.empty() || line[0] == '#'){
//       continue;
//     }
//
//     std::smatch match;
//     if (std::regex_match(line, match, dataRegex))
//     {
//       //
//       std::string elem1 = match[1];
//       std::string elem2 = match[2];
//       std::string elem3 = match[3];
//       std::string paramsStr = match[4];
//
//       //
//       std::istringstream paramsStream(paramsStr);
//       TersoffParams params;
//       paramsStream >> params.m >> params.gamma >> params.lambda3 >> params.c >> params.d
//                    >> params.costheta0 >> params.n >> params.beta >> params.lambda2
//                    >> params.B >> params.R >> params.D >> params.lambda1 >> params.A;
//
//       std::cout << "Read parameters: "
//                 << "m = " << params.m << ", "
//                 << "gamma = " << params.gamma << ", "
//                 << "lambda3 = " << params.lambda3 << ", "
//                 << "c = " << params.c << ", "
//                 << "d = " << params.d << ", "
//                 << "costheta0 = " << params.costheta0 << ", "
//                 << "n = " << params.n << ", "
//                 << "beta = " << params.beta << ", "
//                 << "lambda2 = " << params.lambda2 << ", "
//                 << "B = " << params.B << ", "
//                 << "R = " << params.R << ", "
//                 << "D = " << params.D << ", "
//                 << "lambda1 = " << params.lambda1 << ", "
//                 << "A = " << params.A
//                 << std::endl;
//       //
//       tersoffData[{elem1, elem2, elem3}] = params;
//     }
//   }
//
//   file.close();
// }

// void readTersoffParams(std::ifstream& file, std::vector<TersoffParams>& paramsList)
// {
//
//   std::string line;
//   while (std::getline(file, line)) {
//     //
//     if (line.empty() || line[0] == '#') {
//       continue;
//     }
//
//     //
//     std::istringstream iss(line);
//     double m, gamma, lambda3, c, d, costheta0, n, beta;
//     double lambda2, B, R, D, lambda1, A;
//     std::string element1, element2, element3;
//
//     try
//     {
//       //
//       if (!(iss >> element1 >> element2 >> element3  >>m >> gamma >> lambda3
//         >> c >> d >> costheta0 >> n >> beta>> lambda2 >> B >> R >> D
//         >> lambda1 >> A))
//         {
//           throw std::runtime_error("Invalid format in line: " + line);
//         }
//
//       //
//       paramsList.emplace_back(element1, element2, element3,m, gamma, lambda3,
//         c, d, costheta0, n, beta,lambda2, B, R, D, lambda1, A);
//     }
//     catch (const std::exception& e) {
//       std::cerr << "Error parsing line: " << line << "\n" << e.what() << std::endl;
//     }
//   }
// }

void TerSoff::ReadPotentialFile_fix(std::ifstream& file)
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
              _params = (TersoffParams *) realloc(_params,_maxparam*sizeof(TersoffParams));
              memset(_params + _nparams, 0, DELTA*sizeof(TersoffParams));
           }

            //
            //auto key = std::make_tuple(iname, jname, kname);
            //TersoffParams params;

            iss >> _params[_nparams].m >> _params[_nparams].gamma >> _params[_nparams].lambda3
                >> _params[_nparams].c >> _params[_nparams].d>> _params[_nparams].costheta0 >>
          _params[_nparams].n >>_params[_nparams].beta>> _params[_nparams].lambda2>>
          _params[_nparams].B >> _params[_nparams].R >>_params[_nparams].D >>
          _params[_nparams].lambda1 >> _params[_nparams].A;

          _params[_nparams].m_int = rbmd::Id(_params[_nparams].m);//

            std::cout << "Read parameters: "
            << "m = " << _params[_nparams].m << ", "
            << "gamma = " << _params[_nparams].gamma << ", "
            << "lambda3 = " << _params[_nparams].lambda3 << ", "
            << "c = " << _params[_nparams].c << ", "
            << "d = " << _params[_nparams].d << ", "
            << "costheta0 = " << _params[_nparams].costheta0 << ", "
            << "n = " << _params[_nparams].n << ", "
            << "beta = " << _params[_nparams].beta << ", "
            << "lambda2 = " << _params[_nparams].lambda2 << ", "
            << "B = " << _params[_nparams].B << ", "
            << "R = " << _params[_nparams].R << ", "
            << "D = " << _params[_nparams].D << ", "
            << "lambda1 = " << _params[_nparams].lambda1 << ", "
            << "A = " << _params[_nparams].A
            << std::endl;

            //
            bool unit_convert_flag =true;
            rbmd::Real conversion_factor=1.0;
            if (unit_convert_flag) {
                _params[_nparams].A *= conversion_factor;
                _params[_nparams].B *= conversion_factor;
            }

            //
            //TersoffData[key] = params;
            } catch (const std::exception& e) {
              std::cerr << "Error parsing line: " << line << "\n"
                        << "Error: " << e.what() << std::endl;
              }
        ++_nparams;
    }
    //
    if (_nparams == 0) {
        throw std::runtime_error("No valid parameters found in the potential file.");
    }

  std::cout<< "_nparams: " <<_nparams<<std::endl;
}

void TerSoff::SetupParams()
{
    int i, j, k, m, n;
    //TersoffData params;
    //
  _elem3param = std::vector<std::vector<std::vector<rbmd::Id>>>(
      _nelements,
      std::vector<std::vector<rbmd::Id>>(
          _nelements,
          std::vector<rbmd::Id>(_nelements, -1)
      ));

    //
    for (i = 0; i < _nelements; ++i) {
        for (j = 0; j < _nelements; ++j) {
            for (k = 0; k < _nelements; ++k) {
                n = -1; //
                for (m = 0; m < _nparams; ++m) {
                    if (_params[m].ielement == i &&
                        _params[m].jelement == j &&
                        _params[m].kelement == k) {
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

  for (int i = 0; i < _nparams; ++i) {
    _params[i].cut = _params[i].R + _params[i].D;
    _params[i].cutsq = _params[i].cut * _params[i].cut;

    if (_params[i].n > 0.0) {
      _params[i].c1 = POW(2.0 * _params[i].n * 1.0e-16, -1.0 / _params[i].n);
      _params[i].c2 = POW(2.0 * _params[i].n * 1.0e-8, -1.0 / _params[i].n);
      _params[i].c3 = 1.0 / _params[i].c2;
      _params[i].c4 = 1.0 / _params[i].c1;
    } else {
      _params[i].c1 = _params[i].c2 = _params[i].c3 = _params[i].c4 = 0.0;
    }
  }


  //
  _cutmax = 0.0;
  for (int i = 0; i < _nparams; ++i) {
    if (_params[i].cut > _cutmax) _cutmax = _params[i].cut;
  }


  for (size_t i = 0; i <_nelements; ++i) {
    for (size_t j = 0; j < _nelements; ++j) {
      for (size_t k = 0; k < _nelements; ++k) {
        std::cout << "_elem3param[" << i << "][" << j << "][" << k << "] = "
                  << _elem3param[i][j][k] << std::endl;
      }
    }
  }

}
