#pragma once
#include "../../common/types.h"
#include "force.h"
#include "model/box.h"
#include "../common/erf_table.h"
#include "neighbor_list/include/neighbor_list/neighbor_list.h"
#include "neighbor_list/include/neighbor_list_builder/full_neighbor_list_builder.h"
class TerSoff : public Force
{
public:
  TerSoff();
  virtual ~TerSoff();

  void Init() override;
  void  Execute() override;
  void  ReadPotentialFile(std::ifstream& file);
  void  Element2Type(int narg, char **arg, bool update_setflag);

  struct TersoffParams {
    double m, gamma, lambda3, c, d, costheta0, n, beta;
    double lambda2, B, R, D, lambda1, A;
  };

private:
  std::shared_ptr<BaseNeighborListBuilder> _rbl_neighbor_list_builder;
  std::shared_ptr<BaseNeighborListBuilder> _neighbor_list_builder;
  std::shared_ptr<NeighborList> _rbl_list;
  std::shared_ptr<NeighborList> _list;

  //energy


  //RBL
  std::string _neighbor_type;
  rbmd::Real _cut_off;

  // 用于存储元素对和对应的参数
  typedef std::map<std::tuple<std::string,
  std::string, std::string>, TersoffParams> TersoffData;
  std::ifstream _potential_file;

  rbmd::Id nelements;        // # of unique elements
  char **elements;      // names of unique elements
  rbmd::Id *map;             // mapping from atom types to elements
  int **setflag;      // 0/1 = whether each i,j has been set

};

