#include "nve_ensemble.h"
#include <chrono>

#include "default_position_controller.h"
#include "default_velocity_controller.h"
#include "lj.h"
#include "lj_cut_coul_kspace.h"
#include "cvff.h"
#include "eam.h"
#include "tersoff.h"
#include "energy_stable_scheme_controller.h"
#include "shake_controller.h"
#include "output/include/Logger.hpp"
NVEensemble::NVEensemble() {
  _position_controller = std::make_shared<DefaultPositionController>();
  _velocity_controller = std::make_shared<DefaultVelocityController>();

  // Unified  Force Field Controller
  static const std::unordered_map<std::string, std::function<std::shared_ptr<Force>()>>
  force_map = {
    {"CVFF", [&]() { return std::make_shared<CVFF>(); }},
    {"LJ/CUT", [&]() { return std::make_shared<LJ>(); }},
    {"LJ/CUT/COUL/LONG", [&]() { return std::make_shared<LJCutCoulKspace>(); }},
    {"EAM", [&]() { return std::make_shared<EAM>(); }},
    {"Tersoff", [&]() { return std::make_shared<TerSoff>(); }}
  };

  //force_type
  auto force_type = DataManager::getInstance().getConfigData()->Get<std::string>
  ("type", "hyper_parameters", "force_field");
  if (auto it = force_map.find(force_type); it != force_map.end())
  {
    _force_controller = it->second();
  }
  else {
    Logger::Instance().error("Unsupported force field type: {}", force_type);
  }

  _energy_stable_scheme_controller = std::make_shared<EnergyStableSchemeController>();
}

void NVEensemble::Init() {
  _position_controller->Init();
  _velocity_controller->Init();
  _energy_stable_scheme_controller->Init();

  _force_controller->Init();
  _force_controller->Execute();

  _neighbor_type =
  DataManager::getInstance().getConfigData()->Get<std::string>(
      "type", "hyper_parameters", "neighbor");
}

void NVEensemble::Presolve() {}

void NVEensemble::Solve() {
  auto start = std::chrono::high_resolution_clock::now();

  _velocity_controller->Update();

  _position_controller->Update();

  bool use_shake = false; //TODO: json file
  if (true == use_shake)
  {
    _shake_controller->ShakeA();
  }

  _force_controller->Execute();

  _velocity_controller->Update();

  if (true == use_shake)
  {
    _shake_controller->ShakeB();
  }

  if (_neighbor_type == "RBL")
  {
    _energy_stable_scheme_controller->Update();
  }

  CHECK_RUNTIME(DEVICESYNC());
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration = end - start;

}

void NVEensemble::Postsolve() {}