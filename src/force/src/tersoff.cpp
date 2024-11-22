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

  // 初始化 map
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
  // 打印结果
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

  // 分割字符串
  std::istringstream iss(potential_elements);
  std::vector<std::string> tokens;
  std::string token;
  while (iss >> token) {
    tokens.push_back(token);
  }

  // 设置 narg
  narg = tokens.size();

  // 分配内存给 arg
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

  // 使用 unordered_map 来存储元素与原子类型的映射
  std::unordered_map<std::string, int> element_to_index;

  // 清空之前的数据
  if (_elements) {
    for (i = 0; i < _nelements; i++) {
      delete[] _elements[i];
    }
    delete[] _elements;
  }

  // 初始化数据
  _elements = new char*[ntypes];
  for (i = 0; i < ntypes; i++) {
    _elements[i] = nullptr;
  }

  _nelements = 0;
  _map[0] = -1;  // 保证类型0为无效

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


  // // 遍历每个输入元素名称并映射
  // for (i = 1; i <= narg; i++)
  // {
  //   std::string entry = arg[i - 1];
  //   if (entry == "NULL") {
  //     _map[i] = -1;  // "NULL" 映射为 -1
  //     continue;
  //   }
  //
  //   // 如果元素已经存在，则直接使用现有的索引
  //   auto it = element_to_index.find(entry);
  //   if (it != element_to_index.end()) {
  //     _map[i] = it->second;
  //   } else {
  //     // 新元素，存入 _elements 和 map
  //     _elements[_nelements] = StrDup(entry);
  //     element_to_index[entry] = _nelements;
  //     _map[i] = _nelements;
  //     _nelements++;
  //   }
  // }

  // 更新 setflag 数组
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
//   std::regex dataRegex(R"((\S+)\s+(\S+)\s+(\S+)\s+([\S\s]+))");  // 匹配元素对和参数
//
//   // 读取文件中的每一行
//   while (std::getline(file, line))
//   {
//     // 去掉前后的空格
//     line = line.substr(line.find_first_not_of(" \t"), line.find_last_not_of(" \t") + 1);
//
//     // 跳过注释行
//     if (line.empty() || line[0] == '#'){
//       continue;
//     }
//
//     std::smatch match;
//     if (std::regex_match(line, match, dataRegex))
//     {
//       // 提取元素对和参数
//       std::string elem1 = match[1];
//       std::string elem2 = match[2];
//       std::string elem3 = match[3];
//       std::string paramsStr = match[4];
//
//       // 将参数字符串分割并转换为浮点数
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
//       // 将元素对和参数存入数据结构
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
//     // 跳过空行和注释行
//     if (line.empty() || line[0] == '#') {
//       continue;
//     }
//
//     // 使用 stringstream 解析一行内容
//     std::istringstream iss(line);
//     double m, gamma, lambda3, c, d, costheta0, n, beta;
//     double lambda2, B, R, D, lambda1, A;
//     std::string element1, element2, element3;
//
//     try
//     {
//       // 按指定顺序解析每一列
//       if (!(iss >> element1 >> element2 >> element3  >>m >> gamma >> lambda3
//         >> c >> d >> costheta0 >> n >> beta>> lambda2 >> B >> R >> D
//         >> lambda1 >> A))
//         {
//           throw std::runtime_error("Invalid format in line: " + line);
//         }
//
//       // 创建 TersoffParams 对象并加入列表
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
      // 去掉前后的空格
      line = line.substr(line.find_first_not_of(" \t"), line.find_last_not_of(" \t") + 1);

      // 跳过注释行
      if (line.empty() || line[0] == '#'){
        continue;
      }

      std::istringstream iss(line);
      try {
            // 解析行数据
            std::string iname, jname, kname;
            iss >> iname >> jname >> kname;

            // 确定 ielement, jelement, kelement 是否在元素列表中
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

            // 构建键值并解析参数
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

            // 单位转换（如需要）
            bool unit_convert_flag =true;
            rbmd::Real conversion_factor=1.0;
            if (unit_convert_flag) {
                _params[_nparams].A *= conversion_factor;
                _params[_nparams].B *= conversion_factor;
            }

            // 存储参数到 map 数据结构
            //TersoffData[key] = params;
            } catch (const std::exception& e) {
              std::cerr << "Error parsing line: " << line << "\n"
                        << "Error: " << e.what() << std::endl;
              }
        ++_nparams;
    }
    // 确保至少有一个参数被读取
    if (_nparams == 0) {
        throw std::runtime_error("No valid parameters found in the potential file.");
    }

  std::cout<< "_nparams: " <<_nparams<<std::endl;
}

void TerSoff::SetupParams()
{
    int i, j, k, m, n;
    //TersoffData params;
    // 分配内存给 elem3param，大小为 n_elements x n_elements x n_elements
  _elem3param = std::vector<std::vector<std::vector<rbmd::Id>>>(
      _nelements,
      std::vector<std::vector<rbmd::Id>>(
          _nelements,
          std::vector<rbmd::Id>(_nelements, -1)
      ));

    // 设置 elem3param，确保每个 (i, j, k) 组合在 params 中有唯一匹配
    for (i = 0; i < _nelements; ++i) {
        for (j = 0; j < _nelements; ++j) {
            for (k = 0; k < _nelements; ++k) {
                n = -1; // 初始为未找到
                for (m = 0; m < _nparams; ++m) {
                    if (_params[m].ielement == i &&
                        _params[m].jelement == j &&
                        _params[m].kelement == k) {
                        if (n >= 0) {
                        throw std::runtime_error(
                            std::string("Duplicate entry in potential file for _elements: ") +
                            _elements[i] + " " + _elements[j] + " " + _elements[k]);
                        }
                        n = m; // 记录匹配的参数索引
                    }
                }
                if (n < 0) {
                  throw std::runtime_error(
                    std::string("Missing entry in potential file for _elements: ") +
                    _elements[i] + " " + _elements[j] + " " + _elements[k]);
                }
                _elem3param[i][j][k] = n; // 存储匹配索引
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


  // 设置 cutmax 为所有参数中最大的 cut 值
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
