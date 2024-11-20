#pragma once
#include "pressure_controller.h"
#include "position_controller.h"
#include "velocity_controller.h"

class NoseHooverPressureController : public PressureController {
public:
  NoseHooverPressureController();
  virtual ~NoseHooverPressureController();

  void Init() override;
  void Update() override;
  void InitialIntegrate();
  void FinalIntegrate();

  /**
   * @brief Calculate the current  pressure
   */
  void ComputePressure() override;
  void Couple();
  void ComputeVirial();

  void ComputeTemperature();
  void Computedof();

  void SetUp();
  void NoseHooverChain();
  void ComputeTempTarget();
  void ComputePressTarget();

  void NHCTempIntegrate();
  void NHCPressIntegrate();

  //
  void NHOmegaDot();
  void NH_V_Press();
  void ResetBox();

  void X2Lamda();
  void Lamda2X();

protected:
  std::shared_ptr<PositionController> _position_controller;
  std::shared_ptr<VelocityController> _velocity_controller;

private:
  rbmd::Real _dt;
  rbmd::Real _nktv2p;
  rbmd::Real _kB;
  rbmd::Real _mvv2e;
  rbmd::Real _fmt2v;

  //
  rbmd::Real* _d_temp_contrib;
  rbmd::Real _temp_sum;
  rbmd::Real _temperature;

  // pressure
  rbmd::Real _pressure_start, _pressure_stop, _pressure_damp,_bulkmodulus;//read

  Real3  _p_current, _dilation; //compute
  Real3 _p_start, _p_stop, _p_damp, _p_target;
  rbmd::Real  _scale_factor, _pressure_coupling;

  Real3 _p_freq;
  rbmd::Real _p_freq_max;
  rbmd::Real _p_hydro; // hydrostatic target pressure

  std::vector<rbmd::Id> _p_flag;    //size =6
  rbmd::Id _pdim; // number of barostatted dims

  std::vector<rbmd::Real> _omega, _omega_dot; //size =6
  std::vector<rbmd::Real> _omega_mass;       //size =6
  rbmd::Real _mtk_term1, _mtk_term2; // Martyna-Tobias-Klein corrections
  rbmd::Id _mtk_flag;              // 0 if using Hoover barostat


  //temperature
  rbmd::Real _t_start, _t_stop, _t_damp;   //read
  rbmd::Real  _t_target, _ke_target,_t_freq;//compute
  rbmd::Real _tdof;

  std::vector<rbmd::Real> _eta, _eta_dot; // chain thermostat for particles
  std::vector<rbmd::Real> _eta_dotdot;
  std::vector<rbmd::Real> _eta_mass;
  rbmd::Id _mtchain;              // length of chain
  rbmd::Id _mtchain_default_flag; // 1 = mtchain is default

  std::vector<rbmd::Real> _etap; // chain thermostat for barostat
  std::vector<rbmd::Real> _etap_dot;
  std::vector<rbmd::Real> _etap_dotdot;
  std::vector<rbmd::Real> _etap_mass;
  rbmd::Id _mpchain; // length of chain

  rbmd::Id _nc_tchain, _nc_pchain;
  rbmd::Real _factor_eta;

  //timesteps
  rbmd::Real _dtv, _dtf, _dthalf, _dt4, _dt8, _dto;
  rbmd::Id _eta_mass_flag;
  rbmd::Real _drag, _tdrag_factor; // drag factor on particle thermostat
  rbmd::Real _pdrag_factor;       // drag factor on barostat
};