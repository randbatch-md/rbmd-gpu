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

  void  ReadPotentialElements(const std::string& potential_elements,
                            int& narg, char*** arg);
  void  Element2Type(int narg, char **arg, bool update_setflag);
  void  ReadPotentialFile(std::ifstream& file);
  void ReadPotentialFile_fix(std::ifstream& file);
  void  SetupParams();

  struct TersoffParams
  {
    rbmd::Id ielement, jelement, kelement;
    rbmd::Real m, gamma, lambda3, c, d, costheta0, n, beta;
    rbmd::Real lambda2, B, R, D, lambda1, A;

    rbmd::Real cut, cutsq;
    rbmd::Real c1, c2, c3, c4;
    rbmd::Id m_int;
    std::string iname ,jname ,kname;
  };

  // struct TersoffParams {
  //   double m, gamma, lambda3, c, d, costheta0, n, beta;
  //   double lambda2, B, R, D, lambda1, A;
  //   std::string element1, element2, element3;
  //   rbmd::Real cut, cutsq;
  //   rbmd::Real c1, c2, c3, c4;
  //   rbmd::Id m_int;
  //   std::string iname ,jname ,kname;
  //
  //   TersoffParams(const std::string& e1, const std::string& e2, const std::string& e3,
  //                 double m_, double gamma_, double lambda3_, double c_, double d_,
  //                 double costheta0_, double n_, double beta_, double lambda2_, double B_,
  //                 double R_, double D_, double lambda1_, double A_)
  //       : element1(e1), element2(e2), element3(e3),m(m_), gamma(gamma_),
  //         lambda3(lambda3_), c(c_), d(d_),costheta0(costheta0_), n(n_), beta(beta_),
  //         lambda2(lambda2_), B(B_),R(R_), D(D_), lambda1(lambda1_), A(A_){}
  // };

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

  TersoffParams*  _params;
  rbmd::Id _nelements;        // # of unique elements
  char**    _elements;      // names of unique elements
  std::vector<rbmd::Id> _map;             // mapping from atom types to elements
  //rbmd::Id** _setflag;      // 0/1 = whether each i,j has been set
  std::vector<std::vector<rbmd::Id>> _setflag;



  //
  //std::vector<rbmd::Id> _elem1param;      // mapping from elements to parameters
  //std::vector<std::vector<rbmd::Id>> _elem2param;     // mapping from element pairs to parameters
  std::vector<std::vector<std::vector<rbmd::Id>>> _elem3param;
  //rbmd::Id*** _elem3param;    // mapping from element triplets to parameters
  rbmd::Id _nparams;          // # of stored parameter sets
  rbmd::Id _maxparam;         // max # of parameter sets

  rbmd::Real _cutmax;   //max cutoff for all elements
};

