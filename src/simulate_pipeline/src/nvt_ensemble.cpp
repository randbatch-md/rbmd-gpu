#include "nvt_ensemble.h"

#include <chrono>  //

#include "berendsen_controller.h"
#include "cvff.h"
#include "data_manager.h"
#include "default_position_controller.h"
#include "default_velocity_controller.h"
#include "eam.h"
#include "langevin_controller.h"
#include "lj.h"
#include "lj_cut_coul_kspace.h"
#include "model/md_data.h"
#include "neighbor_list/include/linked_cell/linked_cell_locator.h"
#include "nose_hoover_controller.h"
#include "output/include/Logger.hpp"
#include "rescale_controller.h"
#include "shake_controller.h"
#include "tersoff.h"
NVTensemble::NVTensemble()
{
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

  // unified temperature controller
  _temp_ctrl_type = DataManager::getInstance().getConfigData()->Get
  <std::string>("temp_ctrl_type", "execution");

  if ("RESCALE" == _temp_ctrl_type) {
    _temperature_controller = std::make_shared<RescaleController>();
  }
  else if ("BERENDSEN" == _temp_ctrl_type) {
    _temperature_controller = std::make_shared<BerendsenController>();
  }
  else if ("LANGEVIN" == _temp_ctrl_type) {
    _temperature_controller = std::make_shared<LangevinController>();
  }
  else if ("NOSE_HOOVER" == _temp_ctrl_type) {
    _NoseHoover_controller = std::make_shared<NoseHooverController>();
  }
  else {
    Logger::Instance().error("\033[31m Unsupported temp_ctrl_type: {}\033[0m", _temp_ctrl_type );
    exit(EXIT_FAILURE); //
  }

  // //
  // _NoseHoover_controller = std::make_shared<NoseHooverController>();

  //shake
  _shake_controller = std::make_shared<ShakeController>();
  _integration_type = DataManager::getInstance().getConfigData()->Get
<std::string>("integration_type", "execution");
  std::remove("output_force1.txt");
  std::remove("output_force2.txt");
}

void NVTensemble::Init() {
  _position_controller->Init();
  _velocity_controller->Init();

  _force_controller->Init();
  _force_controller->Execute();
  _shake_controller->Init();

  if (_temperature_controller) {
    _temperature_controller->Init();
  }

  if(_NoseHoover_controller) {
    _NoseHoover_controller->Init();
  }
}

void NVTensemble::Presolve() {}

void NVTensemble::Solve() {
   bool use_shake = DataManager::getInstance().getConfigData()->GetJudge
    <bool>("fix_shake", "hyper_parameters", "extend");

  if ("leapfrog" == _integration_type) {
    if("NOSE_HOOVER" == _temp_ctrl_type)
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
      auto start = std::chrono::high_resolution_clock::now();

      _velocity_controller->Update();
      _position_controller->Update();
      bool use_shake = DataManager::getInstance().getConfigData()->GetJudge<bool>
      ( "fix_shake", "hyper_parameters", "extend");; //TODO: json file
      if (use_shake)
      {
        _shake_controller->ShakeA();
      }

      _force_controller->Execute();

      if ("LANGEVIN"==DataManager::getInstance().getConfigData()->Get<std::string>
        ("temp_ctrl_type", "execution"))
      {
        _temperature_controller->Update();
      }

      _velocity_controller->Update();

      if (use_shake)
      {
        _shake_controller->ShakeB();
      }

      _temperature_controller->ComputeTemperature();

      if ("LANGEVIN" == DataManager::getInstance().getConfigData()->Get<std::string>
        ("temp_ctrl_type", "execution"))
        return;

      _temperature_controller->Update();

      CHECK_RUNTIME(DEVICESYNC());
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<rbmd::Real> duration = end - start;

    }
  }
  else if ("vv" == _integration_type) {
    _position_controller->Update_vv();
    _force_controller->Execute();
    _velocity_controller->Update_vv();

    _temperature_controller->ComputeTemperature();
    _temperature_controller->Update();
  }
  else if ("rkn2" == _integration_type) {
    _position_controller->Update1();
    _force_controller->Execute();
    _velocity_controller->Update1();

    _position_controller->Update2();
    _force_controller->Execute();
    _velocity_controller->Update2();

    _position_controller->Update3();
    _force_controller->Execute();
    // _velocity_controller->Update3();

    _temperature_controller->ComputeTemperature();
    _temperature_controller->Update();
  }
  else if ("prk4" == _integration_type) {
    _velocity_controller->Update1();
    _position_controller->Update1();
    _force_controller->Execute();

    _velocity_controller->Update2();
    _position_controller->Update2();
    _force_controller->Execute();

    _velocity_controller->Update3();
    _position_controller->Update3();
    _force_controller->Execute();

    _velocity_controller->Update4();
    _position_controller->Update4();
    _force_controller->Execute();

    _temperature_controller->ComputeTemperature();
    _temperature_controller->Update();

  }
  else if ("fr4" == _integration_type) {

    _velocity_controller->Update1();
    _position_controller->Update1();
    _force_controller->Execute();

    _velocity_controller->Update2();
    _position_controller->Update2();
    _force_controller->Execute();

    _velocity_controller->Update3();
    _position_controller->Update3();
    _force_controller->Execute();

    _velocity_controller->Update4();
    _position_controller->Update4();
    _force_controller->Execute();

    _temperature_controller->ComputeTemperature();
    _temperature_controller->Update();

  }
  else if ("rkn3a" == _integration_type) {
    _velocity_controller->Update1();
    _position_controller->Update1();
    _force_controller->Execute();

    _velocity_controller->Update2();
    _position_controller->Update2();
    _force_controller->Execute();

    _velocity_controller->Update3();
    _position_controller->Update3();
    _force_controller->Execute();

    _velocity_controller->Update4();
    _position_controller->Update4();
    _force_controller->Execute();

    _temperature_controller->ComputeTemperature();
    _temperature_controller->Update();
  }
  else if ("rkn3b" == _integration_type) {
    _velocity_controller->Update1();
    _position_controller->Update1();
    _force_controller->Execute();

    _velocity_controller->Update2();
    _position_controller->Update2();
    _force_controller->Execute();

    _velocity_controller->Update3();
    _position_controller->Update3();
    _force_controller->Execute();

    _velocity_controller->Update4();
    _position_controller->Update4();
    _force_controller->Execute();

    _temperature_controller->ComputeTemperature();
    _temperature_controller->Update();
  }
  else if ("rkn3c" == _integration_type) {
    _velocity_controller->Update1();
    _position_controller->Update1();
    _force_controller->Execute();

    _velocity_controller->Update2();
    _position_controller->Update2();
    _force_controller->Execute();

    _velocity_controller->Update3();
    _position_controller->Update3();
    _force_controller->Execute();

    _velocity_controller->Update4();
    _position_controller->Update4();
    _force_controller->Execute();

    _temperature_controller->ComputeTemperature();
    _temperature_controller->Update();
  }
  else if ("Beeman" == _integration_type) {
    _position_controller->Update_Beeman();
    _force_controller->Execute();
    _velocity_controller->Update_Beeman();
    _temperature_controller->ComputeTemperature();
    _temperature_controller->Update();
  }
}

void NVTensemble::Postsolve() {}
