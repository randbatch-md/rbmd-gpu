#pragma once
#include "../../common/types.h"

// struct VelocityParams
// {
//   rbmd::Real c1,c2,c3,c4;
// };

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
  void Updatevl() override;
  void Updatebm() override;
 private:
  rbmd::Real _dt;
  rbmd::Real _fmt2v;
  rbmd::Real _c1 = 0.0;
  rbmd::Real _c2 = 0.0;
  rbmd::Real _c3 = 0.0;
  rbmd::Real _c4 = 0.0;
  thrust::device_vector<rbmd::Real> _d_prev_fx;
  thrust::device_vector<rbmd::Real> _d_prev_fy;
  thrust::device_vector<rbmd::Real> _d_prev_fz;
  rbmd::Id _current_step;
  std::string _init_type;
  std::string  _integration_type;
  // rbmd::Real c_type;
  // VelocityParams  _h_v_paras;
  // thrust::device_vector<VelocityParams> _d_v_paras;
};
