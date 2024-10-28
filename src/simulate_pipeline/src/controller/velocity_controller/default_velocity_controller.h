#pragma once
#include "simulate.h"
#include "velocity_controller.h"
class DefaultVelocityController : public VelocityController {
 public:
  DefaultVelocityController();
  virtual ~DefaultVelocityController() = default;

  void Init() override;
  void Update() override;

 private:
  rbmd::Real _dt;
  rbmd::Real _fmt2v;
};
