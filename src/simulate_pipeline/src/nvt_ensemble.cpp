#include "nvt_ensemble.h"

#include <chrono>  // 添加计时功能的库

#include "berendsen_controller.h"
#include "default_position_controller.h"
#include "default_velocity_controller.h"
#include "ljforce.h"
#include "lj_cut_coul_kspace_force.h"
#include "rescale_controller.h"
#include "berendsen_controller.h"
#include "nose_hoover_controller.h"
#include "shake_controller.h"

NVTensemble::NVTensemble()
{
	_position_controller = std::make_shared<DefaultPositionController>(); 
	_velocity_controller = std::make_shared<DefaultVelocityController>(); 
	_force_controller = std::make_shared<LJCutCoulKspaceForce>(); // todo: force_type  json file
	_temperature_controller = std::make_shared<BerendsenController>();
    _shake_controller = std::make_shared<ShakeController>();

}

void NVTensemble::Init() {
  _position_controller->Init();
  _velocity_controller->Init();
  _temperature_controller->Init();
  bool use_shake = DataManager::getInstance().getConfigData()->GetJudge<bool>( "fix_shake", "hyper_parameters", "extend");
  if (use_shake)
  {
      _shake_controller->Init();
  }
  _force_controller->Init();
  _force_controller->Execute();
}

void NVTensemble::Presolve() {}

void NVTensemble::Solve() {
  auto start = std::chrono::high_resolution_clock::now();

  _velocity_controller->Update();

  _position_controller->Update();

  bool use_shake = DataManager::getInstance().getConfigData()->GetJudge<bool>( "fix_shake", "hyper_parameters", "extend");
  if (use_shake)
  {
    _shake_controller->ShakeA();
  }

  _force_controller->Execute();

  _velocity_controller->Update();

  if (use_shake)
  {
    _shake_controller->ShakeB();
  }
  _temperature_controller->Update();

  CHECK_RUNTIME(hipDeviceSynchronize());
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration = end - start;

  std::cout << "time pre step "<< duration.count() << "秒" << std::endl;
}

void NVTensemble::Postsolve() {}
