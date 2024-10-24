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
	//_force_controller = std::make_shared<LJCutCoulKspaceForce>(); // TODO: json file forcetype

    _force_controller = std::make_shared<LJForce>(); // TODO: json file forcetype
	_temperature_controller = std::make_shared<BerendsenController>();
}

void NVTensemble::Init() {
  _position_controller->Init();
  _velocity_controller->Init();
  _temperature_controller->Init();

  _force_controller->Init();
  _force_controller->Execute();
}

void NVTensemble::Presolve() {
  // ����Զ����ʱ Ҫ������Ӧ����
}

void NVTensemble::Solve() {
  auto start = std::chrono::high_resolution_clock::now();

  _velocity_controller->Update();

  _position_controller->Update();

  bool use_shake = false; //TODO: json file 
  if (true == use_shake)
  {
    _shake_controller->ShakeA();
  }

  _force_controller->Execute();

  if ("LANGEVIN"==DataManager::getInstance().getConfigData()->Get<std::string>("temp_ctrl_type", "execution"))
  {
	  _temperature_controller->Update();
  }

  _velocity_controller->Update();

  if (true == use_shake)
  {
    _shake_controller->ShakeB();
  }

  _temperature_controller->ComputeTemp();

  if ("LANGEVIN" == DataManager::getInstance().getConfigData()->Get<std::string>("temp_ctrl_type", "execution"))
	  return;

  _temperature_controller->Update();

  CHECK_RUNTIME(hipDeviceSynchronize());
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration = end - start;

  std::cout << "time pre step "<< duration.count() << "秒" << std::endl;
}

void NVTensemble::Postsolve() {}
