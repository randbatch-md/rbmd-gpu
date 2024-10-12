#include "nvt_ensemble.h"

#include <chrono>  // 添加计时功能的库

#include "berendsen_controller.h"
#include "default_position_controller.h"
#include "default_velocity_controller.h"
#include "ljforce.h"
#include "nose_hoover_controller.h"
#include "rescale_controller.h"
#include "shake_controller.h"

NVTensemble::NVTensemble() {
  _position_controller = std::make_shared<DefaultPositionController>();
  _velocity_controller = std::make_shared<DefaultVelocityController>();
  _force_controller = std::make_shared<LJForce>();  // todo �Զ���forcetype =
  _temperature_controller = std::make_shared<RescaleController>();
}

void NVTensemble::Init() {
  // ���Լ�������ı����ڸ��Ե�init�����ʼ����
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

  bool use_shake;
  if (true == use_shake) {
    _shake_controller->ShakeA();
  }

  _force_controller->Execute();

  _velocity_controller->Update();

  if (true == use_shake) {
    _shake_controller->ShakeB();
  }
  _temperature_controller->Update();

  CHECK_RUNTIME(hipDeviceSynchronize());
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration = end - start;

  std::cout << "time pre step " << duration.count() << "秒" << std::endl;
}

void NVTensemble::Postsolve() {}
