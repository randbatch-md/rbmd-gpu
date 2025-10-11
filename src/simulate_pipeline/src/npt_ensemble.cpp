#include "npt_ensemble.h"

#include <chrono>  //

#include "default_position_controller.h"
#include "default_velocity_controller.h"
#include "lj.h"
#include "lj_cut_coul_kspace.h"
#include "cvff.h"
#include "tersoff.h"
#include "eam.h"
#include "rescale_controller.h"
#include "berendsen_controller.h"
#include "berendsen_pressure_controller.h"
#include "nose_hoover_controller.h"
#include "shake_controller.h"
#include "thermo_stats.hpp"

NPTensemble::NPTensemble() {
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

  // unified temperature/pressure controller
  _temp_ctrl_type = DataManager::getInstance().getConfigData()->Get
    <std::string>("temp_ctrl_type", "execution");
  _press_ctrl_type = DataManager::getInstance().getConfigData()->Get
    <std::string>("press_ctrl_type", "execution");
  if ("BERENDSEN" == _temp_ctrl_type && "BERENDSEN" == _press_ctrl_type) {
    _temperature_controller = std::make_shared<BerendsenController>();
    _pressure_controller = std::make_shared<BerendsenPressureController>();
  }
  else if ("NOSE_HOOVER" == _temp_ctrl_type && "NOSE_HOOVER" == _press_ctrl_type) {
    _NoseHoover_controller = std::make_shared<NoseHooverController>();
  }
  else {
  Logger::Instance().error("\033[31m When ensemble is set to NPT, "
      "temp_ctrl_type must be consistent with press_ctrl_type.\033[0m");
    exit(EXIT_FAILURE); //
  }

  //shake
  _shake_controller = std::make_shared<ShakeController>();
}

void NPTensemble::Init() {
  _position_controller->Init();
  _position_controller->PBC();

  _velocity_controller->Init();

  _force_controller->Init();
  _force_controller->Execute();
  _shake_controller->Init();

  //
  if (_temperature_controller) {
    _temperature_controller->Init();
  }

  if (_pressure_controller) {
    _pressure_controller->Init();
  }
  //
  if (_NoseHoover_controller) {
    _NoseHoover_controller->Init();
  }
}

void NPTensemble::Presolve() {}

void NPTensemble::Solve() {
  auto start = std::chrono::high_resolution_clock::now();

  bool use_shake = DataManager::getInstance().getConfigData()->GetJudge
    <bool>("fix_shake", "hyper_parameters", "extend");

  if("NOSE_HOOVER" == _press_ctrl_type && "NOSE_HOOVER" == _temp_ctrl_type)
  {
    _NoseHoover_controller->InitialIntegrate();//_velocity_controller->Update();
                                           //_position_controller->Update();
    if (true == use_shake)
    {
      _shake_controller->ShakeA();
    }

    _force_controller->Execute();

    _NoseHoover_controller->FinalIntegrate(); //_velocity_controller->Update();

    if (true == use_shake)
    {
      _shake_controller->ShakeB();
    }
  }
  else
  {
     _velocity_controller->Update();     //

     _position_controller->Update();       //

     bool use_shake = false; //TODO: json file
     if (true == use_shake)
     {
       _shake_controller->ShakeA();
     }

     _force_controller->Execute();       //

     _velocity_controller->Update();       //

     if (true == use_shake)
     {
       _shake_controller->ShakeB();
     }

     _temperature_controller->ComputeTemperature();       //

     _temperature_controller->Update();        //
     //
     _pressure_controller->Update();       //
  }

  CHECK_RUNTIME(DEVICESYNC());
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration = end - start;

}

void NPTensemble::Postsolve() {
}