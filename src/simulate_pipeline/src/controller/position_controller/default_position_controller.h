#pragma once
#include <string>

#include "../executioner/executioner.h"
#include "position_controller.h"

class DefaultPositionController : public PositionController {
 public:
  DefaultPositionController();
  virtual ~DefaultPositionController() = default;

  void Init() override;
  void Update() override;

  /**
   * Fill in the center particle coordinates and target particle coordinates of
   * RDF
   */
  void SetCenterTargetPositions();

 private:
  rbmd::Id _num_atoms;
  rbmd::Real _dt;
  std::string _init_type;
};