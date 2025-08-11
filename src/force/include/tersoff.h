#pragma once
#include "../../common/types.h"
#include "force.h"
#include "model/box.h"
#include "../common/erf_table.h"
#include "neighbor_list/include/neighbor_list/neighbor_list.h"
#include "neighbor_list/include/neighbor_list_builder/full_neighbor_list_builder.h"

struct TersoffParams
{
  rbmd::Real m, gamma, lambda3, c, d, costheta0, n, beta;
  rbmd::Real lambda2, B, R, D, lambda1, A;
  rbmd::Real cut, cutsq;
  rbmd::Real c1, c2, c3, c4;
  rbmd::Id ielement, jelement, kelement;
  rbmd::Id m_int;
};

struct ShiftFlag
{
  rbmd::Id   shift_flag;
  rbmd::Real shift_value;
};

class TerSoff : public Force
{
public:
  TerSoff();
  virtual ~TerSoff();

  void Init() override;
  void Execute() override;

  void ReadPotentialElements(const std::string& potential_elements,
                            int& narg, char*** arg);
  void Element2Type(int narg, char **arg, bool update_setflag);
  void ReadPotentialFile(std::ifstream& file);
  void ReadPotentialFile_fix(std::ifstream& file);
  void SetupParams();

  void ComputeTersoff();
  void SumForces();

  void EvaluatePotentialEnergy() override;

private:
  std::shared_ptr<BaseNeighborListBuilder> _rbl_neighbor_list_builder;
  std::shared_ptr<BaseNeighborListBuilder> _neighbor_list_builder;
  std::shared_ptr<NeighborList> _rbl_list;
  std::shared_ptr<NeighborList> _list;
  //energy
  rbmd::Real _e_vdwl = 0;
  rbmd::Real  _e_pe= 0;
  //RBL
  std::string _neighbor_type;
  rbmd::Real _cut_off;

  //
  typedef std::map<std::tuple<std::string,
  std::string, std::string>, TersoffParams> TersoffData;
  std::ifstream _potential_file;

  TersoffParams*  _h_params;
  ShiftFlag _shift;

  rbmd::Id _nelements;        // # of unique elements
  char**    _elements;      // names of unique elements
  std::vector<rbmd::Id> _map;             // mapping from atom types to elements
  //rbmd::Id** _setflag;      // 0/1 = whether each i,j has been set
  std::vector<std::vector<rbmd::Id>> _setflag;

  //
  std::vector<std::vector<std::vector<rbmd::Id>>> _elem3param;

  rbmd::Id _nparams;          // # of stored parameter sets
  rbmd::Id _maxparam;         // max # of parameter sets
  rbmd::Real _cutmax;   //max cutoff for all elements

  //
  thrust::device_vector<rbmd::Id> _d_elem3param;
  thrust::device_vector<rbmd::Id>  _d_map;
  TersoffParams* d_params;
  //
  rbmd::Id  _interval = 2E10;
  std::chrono::duration<rbmd::Real> _duration_list_init;
};

