#pragma once
#include "pressure_controller.h"

class BerendsenPressureController : public PressureController {
public:
  BerendsenPressureController();
  virtual ~BerendsenPressureController();

  void Init() override;
  void Update() override;

  /**
   * @brief Calculate the current  pressure
   */
  void ComputePressure() override;
  void ComputeVirial();

  void Computedof();
  void Couple();
  void ResetBox();
  void X2Lamda();
  void Lamda2X();

private:
  rbmd::Real _dt;
  rbmd::Real _nktv2p;
  rbmd::Real _kB;
  rbmd::Real  _tdof;

  rbmd::Real _pressure_start = 0.0;   //read
  rbmd::Real _pressure_stop  = 0.0;   //read
  rbmd::Real _pressure_damp  = 1000;   //read
  rbmd::Real  _bulkmodulus = 100;       //read

  Real3  _p_current, _dilation;
  Real3 _p_start, _p_stop, _p_damp, _p_target;

};