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

TerSoff::TerSoff()
{
  _rbl_neighbor_list_builder = std::make_shared<RblFullNeighborListBuilder>();
  _neighbor_list_builder = std::make_shared<FullNeighborListBuilder>();
  params = nullptr;
}

TerSoff::~TerSoff()
{
  CHECK_RUNTIME(FREE(params));
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

  //SetupParams();
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
  if (elements) {
    for (i = 0; i < nelements; i++) {
      delete[] elements[i];
    }
    delete[] elements;
  }

  // 初始化数据
  elements = new char*[ntypes];
  for (i = 0; i < ntypes; i++) {
    elements[i] = nullptr;
  }

  nelements = 0;
  map[0] = -1;  // 保证类型0为无效

  // 遍历每个输入元素名称并映射
  for (i = 1; i <= narg; i++)
  {
    std::string entry = arg[i - 1];
    if (entry == "NULL") {
      map[i] = -1;  // "NULL" 映射为 -1
      continue;
    }

    // 如果元素已经存在，则直接使用现有的索引
    auto it = element_to_index.find(entry);
    if (it != element_to_index.end()) {
      map[i] = it->second;
    } else {
      // 新元素，存入 elements 和 map
      elements[nelements] = StrDup(entry);
      element_to_index[entry] = nelements;
      map[i] = nelements;
      nelements++;
    }
  }

  // // 更新 setflag 数组
  // if (update_setflag) {
  //   int count = 0;
  //   for (i = 1; i <= ntypes; i++) {
  //     for (j = i; j <= ntypes; j++) {
  //       setflag[i][j] = 0;
  //       if ((map[i] >= 0) && (map[j] >= 0)) {
  //         setflag[i][j] = 1;
  //         count++;
  //       }
  //     }
  //   }
  //
  //   if (count == 0) {
  //     std::cerr << "Incorrect args for tersoff coefficients"<< std::endl;
  //   }
  // }
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

void TerSoff::ReadPotentialFile_fix(std::ifstream& file)
{
    // 清理现有数据并初始化
    //TersoffData.clear();
    //params = nullptr;
    nparams = maxparam = 0;

    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);

        try {
            // 解析行数据
            std::string iname, jname, kname;
            iss >> iname >> jname >> kname;

            // 确定 ielement, jelement, kelement 是否在元素列表中
            int ielement, jelement, kelement;
            for (ielement = 0; ielement < nelements; ++ielement)
                if (iname == elements[ielement]) break;
            if (ielement == nelements) continue;

            for (jelement = 0; jelement < nelements; ++jelement)
                if (jname == elements[jelement]) break;
            if (jelement == nelements) continue;

            for (kelement = 0; kelement < nelements; ++kelement)
                if (kname == elements[kelement]) break;
            if (kelement == nelements) continue;

            int DELTA =4;
            if (nparams == maxparam)
            {
              maxparam += DELTA;
              params = (TersoffParams *) realloc(params,maxparam*sizeof(TersoffParams));
              memset(params + nparams, 0, DELTA*sizeof(TersoffParams));
           }

            // 构建键值并解析参数
            //auto key = std::make_tuple(iname, jname, kname);
            //TersoffParams params;

            iss >> params[nparams].m >> params[nparams].gamma >> params[nparams].lambda3
                >> params[nparams].c >> params[nparams].d>> params[nparams].costheta0 >>
          params[nparams].n >> params[nparams].beta>> params[nparams].lambda2>>
          params[nparams].B >> params[nparams].R >>params[nparams].D >>
          params[nparams].lambda1 >> params[nparams].A;

          std::cout << "Read parameters: "
          << "m = " << params[nparams].m << ", "
          << "gamma = " << params[nparams].gamma << ", "
          << "lambda3 = " << params[nparams].lambda3 << ", "
          << "c = " << params[nparams].c << ", "
          << "d = " << params[nparams].d << ", "
          << "costheta0 = " << params[nparams].costheta0 << ", "
          << "n = " << params[nparams].n << ", "
          << "beta = " << params[nparams].beta << ", "
          << "lambda2 = " << params[nparams].lambda2 << ", "
          << "B = " << params[nparams].B << ", "
          << "R = " << params[nparams].R << ", "
          << "D = " << params[nparams].D << ", "
          << "lambda1 = " << params[nparams].lambda1 << ", "
          << "A = " << params[nparams].A
          << std::endl;


            // 单位转换（如需要）
            bool unit_convert_flag =true;
            rbmd::Real conversion_factor=1.0;
            if (unit_convert_flag) {
                params[nparams].A *= conversion_factor;
                params[nparams].B *= conversion_factor;
            }

            // 存储参数到 map 数据结构
            //TersoffData[key] = params;
            ++nparams;

        } catch (const std::exception& e) {
            std::cerr << "Error parsing line: " << line << "\n"
                      << "Error: " << e.what() << std::endl;
            continue; // 跳过错误的行
        }
    }

    // 确保至少有一个参数被读取
    if (nparams == 0) {
        throw std::runtime_error("No valid parameters found in the potential file.");
    }

  std::cout<< "nparams: " <<nparams<<std::endl;
}



void TerSoff::SetupParams()
{
    int i, j, k, m, n;
    //TersoffData params;
    // 分配内存给 elem3param，大小为 nelements x nelements x nelements
  //elem3param[nelements][nelements][nelements];

    // 设置 elem3param，确保每个 (i, j, k) 组合在 params 中有唯一匹配
    for (i = 0; i < nelements; ++i) {
        for (j = 0; j < nelements; ++j) {
            for (k = 0; k < nelements; ++k) {
                n = -1; // 初始为未找到
                for (m = 0; m < nparams; ++m) {
                    if (params[m].ielement == i &&
                        params[m].jelement == j &&
                        params[m].kelement == k) {
                        if (n >= 0) {
                        throw std::runtime_error(
                            std::string("Duplicate entry in potential file for elements: ") +
                            elements[i] + " " + elements[j] + " " + elements[k]);
                        }
                        n = m; // 记录匹配的参数索引
                    }
                }
                if (n < 0) {
                  throw std::runtime_error(
                    std::string("Missing entry in potential file for elements: ") +
                    elements[i] + " " + elements[j] + " " + elements[k]);
                }
                elem3param[i][j][k] = n; // 存储匹配索引
            }
        }
    }

    // 计算 derived 参数，例如 cut 和 cutsq
  for (int i = 0; i < nparams; ++i) {
    params[i].cut = params[i].R + params[i].D;
    params[i].cutsq = params[i].cut * params[i].cut;

    if (params[i].n > 0.0) {
      params[i].c1 = POW(2.0 * params[i].n * 1.0e-16, -1.0 / params[i].n);
      params[i].c2 = POW(2.0 * params[i].n * 1.0e-8, -1.0 / params[i].n);
      params[i].c3 = 1.0 / params[i].c2;
      params[i].c4 = 1.0 / params[i].c1;
    } else {
      params[i].c1 = params[i].c2 = params[i].c3 = params[i].c4 = 0.0;
    }
  }


  // 设置 cutmax 为所有参数中最大的 cut 值
  cutmax = 0.0;
  for (int i = 0; i < nparams; ++i) {
    if (params[i].cut > cutmax) cutmax = params[i].cut;
  }


  for (size_t i = 0; i <nelements; ++i) {
    for (size_t j = 0; j < nelements; ++j) {
      for (size_t k = 0; k < nelements; ++k) {
        std::cout << "elem3param[" << i << "][" << j << "][" << k << "] = "
                  << elem3param[i][j][k] << std::endl;
      }
    }
  }

}
