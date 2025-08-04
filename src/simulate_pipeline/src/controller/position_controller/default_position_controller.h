#pragma once
#include <string>

#include "position_controller.h"

class DefaultPositionController : public PositionController {
 public:
  DefaultPositionController();
  virtual ~DefaultPositionController() = default;

  void Init() override;
  void Update() override;
  void Update1() override;
  void Update2() override;
  void Update3() override;
  void Update4() override;
  void Update_vv() override;


  /**
   * Fill in the center particle coordinates and target particle coordinates of
   * RDF
   */
  void SetCenterTargetPositions();

 private:
  rbmd::Real _dt;
  std::string _init_type;
  rbmd::Real _fmt2v;
  std::string  _integration_type;
  rbmd::Real _d1 = 0.0;
  rbmd::Real _d2 = 0.0;
  rbmd::Real _d3 = 0.0;
  rbmd::Real _d4 = 0.0;
};
