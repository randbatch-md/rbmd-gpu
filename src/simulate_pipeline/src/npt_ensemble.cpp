#include "npt_ensemble.h"

#include <chrono>  // 添加计时功能的库

#include "default_position_controller.h"
#include "default_velocity_controller.h"
#include "lj.h"
#include "lj_cut_coul_kspace.h"
#include "cvff.h"
#include "tersoff.h"
#include "rescale_controller.h"
#include "berendsen_controller.h"
#include "berendsen_pressure_controller.h"
#include "nose_hoover_controller.h"
#include "shake_controller.h"

NPTensemble::NPTensemble() {
  //_position_controller = std::make_shared<DefaultPositionController>();
  //_velocity_controller = std::make_shared<DefaultVelocityController>();
  _force_controller = std::make_shared<LJCutCoulKspace>(); // TODO: json file forcetype
  //_temperature_controller = std::make_shared<BerendsenController>();
  //_pressure_controller = std::make_shared<BerendsenPressureController>();
  _NoseHoover_controller = std::make_shared<NoseHooverPressureController>();

  _shake_controller = std::make_shared<ShakeController>();
}

void NPTensemble::Init() {
  //_position_controller->Init();
  //_velocity_controller->Init();
  //_temperature_controller->Init();
  //_pressure_controller->Init();
  _NoseHoover_controller->Init();

  _force_controller->Init();
  _force_controller->Execute();

  _shake_controller->Init();
}

void NPTensemble::Presolve() {}

void NPTensemble::Solve() {
  auto start = std::chrono::high_resolution_clock::now();

  bool use_shake = DataManager::getInstance().getConfigData()->GetJudge
    <bool>("fix_shake", "hyper_parameters", "extend");
  auto press_ctrl_type = DataManager::getInstance().getConfigData()->Get
    <std::string>("press_ctrl_type", "execution");
  auto temp_ctrl_type = DataManager::getInstance().getConfigData()->Get
  <std::string>("temp_ctrl_type", "execution");

  if("NOSE_HOOVER" == temp_ctrl_type)
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
  // else if("BERENDSEN" == press_ctrl_type)
  // {
  //   _velocity_controller->Update();
  //
  //   _position_controller->Update();
  //
  //   bool use_shake = false; //TODO: json file
  //   if (true == use_shake)
  //   {
  //     _shake_controller->ShakeA();
  //   }
  //
  //   _force_controller->Execute();
  //
  //   if ("LANGEVIN"==DataManager::getInstance().getConfigData()->
  //     Get<std::string>("temp_ctrl_type", "execution"))
  //   {
  //     _temperature_controller->Update();
  //   }
  //
  //   _velocity_controller->Update();
  //
  //   if (true == use_shake)
  //   {
  //     _shake_controller->ShakeB();
  //   }
  //
  //   _temperature_controller->ComputeTemperature();
  //
  //   if ("LANGEVIN" == DataManager::getInstance().getConfigData()->Get<std::string>("temp_ctrl_type", "execution"))
  //     return;
  //
  //   _temperature_controller->Update();
  //   //
  //   _pressure_controller->Update();
  //}


  CHECK_RUNTIME(hipDeviceSynchronize());
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration = end - start;

  std::cout << "time pre step "<< duration.count() << "秒" << std::endl;
}

void NPTensemble::Postsolve() {
}