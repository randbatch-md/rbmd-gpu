#include "nvt_ensemble.h"

#include <chrono>  //

#include "data_manager.h"
#include "model/md_data.h"

#include "default_position_controller.h"
#include "default_velocity_controller.h"
#include "berendsen_controller.h"
#include "langevin_controller.h"
#include "rescale_controller.h"
#include "nose_hoover_controller.h"
#include "cvff.h"
#include "lj_cut_coul_kspace.h"
#include "lj.h"
#include "tersoff.h"
#include "eam.h"

#include "shake_controller.h"
#include "output/include/Logger.hpp"
#if USE_MACE
   #include "maceload.h"
#endif

NVTensemble::NVTensemble()
{
  _position_controller = std::make_shared<DefaultPositionController>();
  _velocity_controller = std::make_shared<DefaultVelocityController>();

  // Unified  Force Field Controller
  #if USE_MACE
     static const std::unordered_map<std::string, std::function<std::shared_ptr<Force>()>>
     force_map = {
       {"CVFF", [&]() { return std::make_shared<CVFF>(); }},
       {"LJ/CUT", [&]() { return std::make_shared<LJ>(); }},
       {"LJ/CUT/COUL/LONG", [&]() { return std::make_shared<LJCutCoulKspace>(); }},
       {"EAM", [&]() { return std::make_shared<EAM>(); }},
       {"Tersoff", [&]() { return std::make_shared<TerSoff>(); }},
       {"MACE", [&]() { return std::make_shared<maceload>(); }}
     };
  //    std::cout<<"mace"<<std::endl;
  #else
     static const std::unordered_map<std::string, std::function<std::shared_ptr<Force>()>>
     force_map = {
       {"CVFF", [&]() { return std::make_shared<CVFF>(); }},
       {"LJ/CUT", [&]() { return std::make_shared<LJ>(); }},
       {"LJ/CUT/COUL/LONG", [&]() { return std::make_shared<LJCutCoulKspace>(); }},
       {"EAM", [&]() { return std::make_shared<EAM>(); }},
       {"Tersoff", [&]() { return std::make_shared<TerSoff>(); }}
     };
  //     std::cout<<"nomace"<<std::endl;
  #endif
  //force_type
  auto force_type = DataManager::getInstance().getConfigData()->Get<std::string>
  ("type", "hyper_parameters", "force_field");
  std::cout<<force_type<<std::endl;
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

void NVTensemble::Postsolve() {}
