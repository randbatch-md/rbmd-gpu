#include "nvt_ensemble.h"

#include <chrono>  // 添加计时功能的库

#include "default_position_controller.h"
#include "default_velocity_controller.h"
#include "lj.h"
#include "lj_cut_coul_kspace.h"
#include "cvff.h"
#include "tersoff.h"
#include "rescale_controller.h"
#include "berendsen_controller.h"
#include "nose_hoover_controller.h"
#include "shake_controller.h"
#include "data_manager.h"
#include "model/md_data.h"
#include "maceload.h"

NVTensemble::NVTensemble()
{
  _position_controller = std::make_shared<DefaultPositionController>();
  _velocity_controller = std::make_shared<DefaultVelocityController>();

  auto force_type = DataManager::getInstance().getConfigData()->Get<std::string>
    ("type", "hyper_parameters", "force_field");
  if ("CVFF" == force_type) {
    _force_controller = std::make_shared<CVFF>(); // TODO: json file forcetype
  }
  else if ("LJ/CUT" == force_type){
    _force_controller = std::make_shared<LJ>(); // TODO: json file forcetype
  }
  else if ("LJ/CUT/COUL/LONG" == force_type){
    _force_controller = std::make_shared<LJCutCoulKspace>(); // TODO: json file forcetype
  }
  else if ("MACE" == force_type) {
    std::cout<<"mace_nvt_run"<<std::endl;
    _force_controller = std::make_shared<maceload>();
  }
  _temperature_controller = std::make_shared<BerendsenController>();
  _shake_controller = std::make_shared<ShakeController>();

  _NoseHoover_controller = std::make_shared<NoseHooverController>();
}

void NVTensemble::Init() {
  _position_controller->Init();
  _velocity_controller->Init();
  _temperature_controller->Init();

  _force_controller->Init();
  _force_controller->Execute();
  _shake_controller->Init();

  _temp_ctrl_type = DataManager::getInstance().getConfigData()->Get
  <std::string>("temp_ctrl_type", "execution");
  if("NOSE_HOOVER" == _temp_ctrl_type) {
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
