#pragma once
#include <string>

#include "position_controller.h"

class DefaultPositionController : public PositionController {
 public:
  DefaultPositionController();
  virtual ~DefaultPositionController() = default;

  void Init() override;
  void Update() override;
  void Updatebm() override;
  /**
   * Fill in the center particle coordinates and target particle coordinates of
   * RDF
   */
  void SetCenterTargetPositions();

 private:
  rbmd::Real _dt;
  std::string _init_type;
  rbmd::Real _fmt2v;

  rbmd::Real _par_a = 0.0;
  rbmd::Real _par_b = 0.0;
  rbmd::Id _current_step;
  thrust::device_vector<rbmd::Real> _d_prev_fx;
  thrust::device_vector<rbmd::Real> _d_prev_fy;
  thrust::device_vector<rbmd::Real> _d_prev_fz;
  thrust::device_vector<rbmd::Real> _d_prev_px;
  thrust::device_vector<rbmd::Real> _d_prev_py;
  thrust::device_vector<rbmd::Real> _d_prev_pz;
};
