#pragma once

#include "velocity_controller.h"
class DefaultVelocityController : public VelocityController {
 public:
  DefaultVelocityController();
  virtual ~DefaultVelocityController() = default;

  void Init() override;
  void Update() override;
  void Update1() override;
  void Update2() override;
  void Update3() override;
  void Update4() override;
  void Update_vv() override;
  void Update_Beeman() override;


 private:
  rbmd::Real _dt;
  rbmd::Real _fmt2v;
  std::string  _integration_type;
  rbmd::Real _c1 = 0.0;
  rbmd::Real _c2 = 0.0;
  rbmd::Real _c3 = 0.0;
  rbmd::Real _c4 = 0.0;
  thrust::device_vector<rbmd::Real> _d_prev_fx;
  thrust::device_vector<rbmd::Real> _d_prev_fy;
  thrust::device_vector<rbmd::Real> _d_prev_fz;
  thrust::device_vector<rbmd::Real>_d_prev_px;
  thrust::device_vector<rbmd::Real>_d_prev_py;
  thrust::device_vector<rbmd::Real>_d_prev_pz;
};
