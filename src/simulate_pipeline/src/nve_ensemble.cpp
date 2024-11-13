#include "nve_ensemble.h"
#include <chrono>

#include "default_position_controller.h"
#include "default_velocity_controller.h"
#include "ljforce.h"
#include "lj_cut_coul_kspace_force.h"
#include "cvff.h"
#include "tersoff.h"
#include "energy_stable_scheme_controller.h"
#include "shake_controller.h"

NVEensemble::NVEensemble() {
  _position_controller = std::make_shared<DefaultPositionController>();
  _velocity_controller = std::make_shared<DefaultVelocityController>();
  _force_controller = std::make_shared<LJCutCoulKspaceForce>(); // TODO: json file forcetype
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

  CHECK_RUNTIME(hipDeviceSynchronize());
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration = end - start;

  std::cout << "time pre step "<< duration.count() << "秒" << std::endl;
}

void NVEensemble::Postsolve() {}