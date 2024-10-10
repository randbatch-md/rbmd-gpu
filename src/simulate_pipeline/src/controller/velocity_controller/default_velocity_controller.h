#pragma once
#include "../executioner/executioner.h"
#include "velocity_controller.h"
class DefaultVelocityController : public VelocityController {
 public:
  DefaultVelocityController();
  virtual ~DefaultVelocityController() = default;

  void Init() override;
  void Update() override;

 private:
  rbmd::Id _num_atoms;
  rbmd::Real _dt;
  rbmd::Real _fmt2v;
};