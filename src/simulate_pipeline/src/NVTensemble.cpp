#include "NVTensemble.h"
#include "controller/position_controller/default_position_controller.h"
#include "controller/velocity_controller/default_velocity_controller.h"
#include "controller/shake_controller/shake_controller.h"
#include "controller/force_controller/default_force_controller.h"
#include "controller/temperature_controller/rescale_controller.h"

NVTensemble::NVTensemble()
{
	_position_controller = std::make_shared<DefaultPositionController>(); 
	_velocity_controller = std::make_shared<DefaultVelocityController>(); 
	_force_controller = std::make_shared<DefaultForceController>(); 
	_temperature_controller = std::make_shared<RescaleController>();

}

void NVTensemble::Init()
{
	// 各自计算所需的变量在各自的init里面初始化；
	_position_controller->Init();
	_velocity_controller->Init();
	_temperature_controller->Init();

	_force_controller->Init(); 
	_force_controller->Execute();
}

void NVTensemble::Presolve()
{
	// 计算远程力时 要添加相应参数
}

void NVTensemble::Solve()
{
	_velocity_controller->Update();
	
	_position_controller->Update();
	
	bool use_shake;
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
	_temperature_controller->Update();
}

void NVTensemble::Postsolve()
{

}
