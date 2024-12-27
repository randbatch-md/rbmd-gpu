#pragma once
#include "../../common/types.h"
#include "force.h"
#include "model/box.h"
#include "../common/erf_table.h"
#include "neighbor_list/include/neighbor_list/neighbor_list.h"
#include "neighbor_list/include/neighbor_list_builder/full_neighbor_list_builder.h"
class CVFF : public Force
{
public:
  CVFF();
  virtual ~CVFF();

  void Init() override;
  void  Execute() override;
  void EvaluatePotentialenergy() override;

  void ComputeLJCutCoulForce();
  void ComputeLJVerlet() ;
  void ComputeLJRBL();

  void ComputeKspaceForce();
  void SumForces();

  //
  void ComputeQsqSum();
  rbmd::Real ComputeRMS(rbmd::Id kmax,rbmd::Real box_length,rbmd::Real q2);
  void SetKspacePara();

  void  ComputeChargeStructureFactorEwald(
          Box box,
          rbmd::Id num_atoms,
          Int3 Kmax_array,
          rbmd::Real alpha,
          rbmd::Real qqr2e,
          rbmd::Real* value_Re_array,
          rbmd::Real* value_Im_array);
  void ComputeEwlad();//Ewald

  void RBEInit(Box box,rbmd::Real alpha,rbmd::Id RBE_P);
  void GetPsampleKey();
  void ComputeChargeStructureFactorRBE(
         Box box,
         rbmd::Id num_atoms,
         Int3 Kmax_array,
         rbmd::Real alpha,
         rbmd::Id RBE_P,
         rbmd::Real qqr2e,
         thrust::device_vector<rbmd::Real> rhok_real_redue,
         thrust::device_vector<rbmd::Real> rhok_image_redue);
  void ComputeRBE();//RBE


  void ComputeLJCoulEnergy();

  void ComputeSelfEnergy(
    rbmd::Real alpha,
    rbmd::Real qqr2e,
    rbmd::Real& ave_self_energy);  //self  Energy

  void ComputeKspaceEnergy(
        Box box,
        rbmd::Id _num_atoms,
        Int3 Kmax_array,
        rbmd::Real alpha,
        rbmd::Real qqr2e,
        rbmd::Real& ave_ekspace);   //Ewald  Energy

  void ComputeBondForce(); //Harmonic
  void ComputeAngleForce(); //Harmonic
  void ComputeDihedralForce(); //Harmonic
  void ComputeImproperForce(); //Harmonic

private:
  std::shared_ptr<BaseNeighborListBuilder> _rbl_neighbor_list_builder;
  std::shared_ptr<BaseNeighborListBuilder> _neighbor_list_builder;
  std::shared_ptr<NeighborList> _rbl_list;
  std::shared_ptr<NeighborList> _list;

  //energy
  rbmd::Real _ave_evdwl = 0;
  rbmd::Real _ave_ecoul = 0;
  rbmd::Real _ave_especial_coul = 0;
  rbmd::Real _ave_self_energy = 0;
  rbmd::Real _ave_ekspace = 0;
  rbmd::Real _ave_ebond = 0;
  rbmd::Real _ave_eangle = 0;
  rbmd::Real _ave_edihedral = 0;
  rbmd::Real _ave_pe = 0;

  rbmd::Real _ave_evdwl_rbl = 0;
  rbmd::Real _ave_ecoul_rbl = 0;
  rbmd::Real _ave_pe_rbl = 0;

  //RBL
  std::string _neighbor_type;
  rbmd::Real _cut_off;

  rbmd::Real _corr_value_x = 0;
  rbmd::Real _corr_value_y = 0;
  rbmd::Real _corr_value_z = 0;

  //EWALD
  rbmd::Real _accuracy;
  rbmd::Real _g_ewald;
  rbmd::Real _alpha;
  rbmd::Real _sum_sq_charge;
  rbmd::Real* _h_Re_array;
  rbmd::Real* _h_Im_array;

  rbmd::Id _Kmax;
  rbmd::Id _Kmax3D;
  rbmd::Id kmax_x, kmax_y, kmax_z;
  Int3 _kmax_array;
  Real3 _unitk;
  rbmd::Real _gsqmx;
  rbmd::Id _kmax_x_orig, _kmax_y_orig, _kmax_z_orig;
  rbmd::Id kcount;

  //RBE
  std::string _coulomb_type;
  rbmd::Id _RBE_P;
  rbmd::Real _qqr2e;
  rbmd::Id _num_k;
  thrust::device_vector<rbmd::Real>  _P_Sample_x;
  thrust::device_vector<rbmd::Real>  _P_Sample_y;
  thrust::device_vector<rbmd::Real>  _P_Sample_z;
  //
  thrust::device_vector<rbmd::Id>  _psample_key;
  thrust::device_vector<rbmd::Real> _rhok_real_redue;
  thrust::device_vector<rbmd::Real> _rhok_image_redue;

};

