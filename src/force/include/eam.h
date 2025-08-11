#pragma once
#include "../../common/types.h"
#include "force.h"
#include "model/box.h"
#include "neighbor_list/include/neighbor_list/neighbor_list.h"
#include "neighbor_list/include/neighbor_list_builder/full_neighbor_list_builder.h"

struct EAMParameters {

  rbmd::Real drho, dr, rhomax, rhomin;;
  rbmd::Id nrho, nr;
};

class EAM : public Force
{
public:
  EAM();
  virtual ~EAM();

  void Init() override;
  void  Execute() override;
  void EvaluatePotentialEnergy() override;

  void EAMVerlet();
  void EAMRBL();

  void SumForces();
  void ReadPotentialFile(const std::string& filename);
  void InitStyle();
  void AllocateEAM();
  void file2array();
  void interpolate(rbmd::Id n, rbmd::Real delta, std::vector<rbmd::Real>& f,  thrust::host_vector<Real7>& spline);
  void array2spline();
  void SetEAM();

  void ComputEAMEnergy();

private:
  std::shared_ptr<BaseNeighborListBuilder> _rbl_neighbor_list_builder;
  std::shared_ptr<BaseNeighborListBuilder> _neighbor_list_builder;
  std::shared_ptr<NeighborList> _rbl_list;
  std::shared_ptr<NeighborList> _list;

  //energy

  rbmd::Real _e_embedding = 0;
  rbmd::Real _e_pair = 0;
  rbmd::Real _e_pe = 0;

  //RBL
  std::string  _energy_rbl_flag = "yes";
  std::string _neighbor_type;
  rbmd::Real _cut_off;

  rbmd::Real _corr_value_x = 0;
  rbmd::Real _corr_value_y = 0;
  rbmd::Real _corr_value_z = 0;


  struct Funcfl
  {
    rbmd::Id nrho, nr;
    rbmd::Real drho, dr, cut_off;
    std::vector<rbmd::Real> frho;
    std::vector<rbmd::Real> zr;
    std::vector<rbmd::Real> rhor;
  };
  Funcfl file;

  EAMParameters eam_paras;

  std::vector<rbmd::Real> frho;
  std::vector<rbmd::Real> z2r;
  std::vector<rbmd::Real> rhor;

  std::vector<rbmd::Id> type2frho;
  std::vector<Id2> type2rhor;
  std::vector<Id2> type2z2r;
  std::vector<Real2> scale;


  thrust::host_vector<Real7> _h_frho_spline;
  thrust::host_vector<Real7> _h_rhor_spline;
  thrust::host_vector<Real7> _h_z2r_spline;

  thrust::device_vector<Real7> _d_frho_spline;
  thrust::device_vector<Real7> _d_rhor_spline;
  thrust::device_vector<Real7> _d_z2r_spline;

};

