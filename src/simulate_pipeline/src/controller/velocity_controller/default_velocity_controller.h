#pragma once

#include "velocity_controller.h"
class DefaultVelocityController : public VelocityController {
 public:
  DefaultVelocityController();
  virtual ~DefaultVelocityController() = default;

  void Init() override;
  void Update() override;
  void Updatebm() override;
 private:
  rbmd::Real _dt;
  rbmd::Real _fmt2v;

  rbmd::Real _par_a = 1.0;
  rbmd::Real _par_b = 3.0;
  thrust::device_vector<rbmd::Real> _d_prev_fx;
  thrust::device_vector<rbmd::Real> _d_prev_fy;
  thrust::device_vector<rbmd::Real> _d_prev_fz;
  thrust::device_vector<rbmd::Real> _d_pr_prev_fx;
  thrust::device_vector<rbmd::Real> _d_pr_prev_fy;
  thrust::device_vector<rbmd::Real> _d_pr_prev_fz;
};
