#pragma once
#include "types.h"
#include "vector"
class ShakeController {
 public:
  ShakeController();
  virtual ~ShakeController() = default;

  void ShakeA();
  void ShakeB();

 public:
  std::vector<rbmd::Real> _shake_px;
  std::vector<rbmd::Real> _shake_py;
  std::vector<rbmd::Real> _shake_pz;

  std::vector<rbmd::Real> _shake_vx;
  std::vector<rbmd::Real> _shake_vy;
  std::vector<rbmd::Real> _shake_vz;
};