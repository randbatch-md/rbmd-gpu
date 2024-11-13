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
}

TerSoff::~TerSoff()
{
}

void TerSoff::Init()
{
  auto potential_file = DataManager::getInstance().getConfigData()->Get
 <std::string>("potential_file", "hyper_parameters", "force_field");
  _potential_file.open(potential_file);
  ReadPotentialFile(_potential_file);
}

void TerSoff::Execute()
{
}

void TerSoff::ReadPotentialFile(std::ifstream& file)
{
  if (!file.is_open())
  {
    std::cerr << "Unable to open the file." << std::endl;
  }

  TersoffData tersoffData;
  std::string line;
  std::regex dataRegex(R"((\S+)\s+(\S+)\s+(\S+)\s+([\S\s]+))");  // 匹配元素对和参数

  // 读取文件中的每一行
  while (std::getline(file, line))
  {
    // 去掉前后的空格
    line = line.substr(line.find_first_not_of(" \t"), line.find_last_not_of(" \t") + 1);

    // 跳过注释行
    if (line.empty() || line[0] == '#'){
      continue;
    }

    std::smatch match;
    if (std::regex_match(line, match, dataRegex))
    {
      // 提取元素对和参数
      std::string elem1 = match[1];
      std::string elem2 = match[2];
      std::string elem3 = match[3];
      std::string paramsStr = match[4];

      // 将参数字符串分割并转换为浮点数
      std::istringstream paramsStream(paramsStr);
      TersoffParams params;
      paramsStream >> params.m >> params.gamma >> params.lambda3 >> params.c >> params.d
                   >> params.costheta0 >> params.n >> params.beta >> params.lambda2
                   >> params.B >> params.R >> params.D >> params.lambda1 >> params.A;

      // 将元素对和参数存入数据结构
      tersoffData[{elem1, elem2, elem3}] = params;
    }
  }

  file.close();
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

  // 更新 setflag 数组
  if (update_setflag) {
    int count = 0;
    for (i = 1; i <= ntypes; i++) {
      for (j = i; j <= ntypes; j++) {
        setflag[i][j] = 0;
        if ((map[i] >= 0) && (map[j] >= 0)) {
          setflag[i][j] = 1;
          count++;
        }
      }
    }

    if (count == 0) {
      std::cerr << "Incorrect args for tersoff coefficients"<< std::endl;
    }
  }
}
