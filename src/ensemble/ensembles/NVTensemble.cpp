#include "NVTensemble.h"
#include "default_position_controller.h"
#include "default_velocity_controller.h"

NVTensemble::NVTensemble()
{
	_position_controller=std::make_shared<DefaultPositionController>(); 
	_velocity_controller=std::make_shared<DefaultVelocityController>(); 
	_force_controller = std::make_shared<DefaultForceController>();
}

void NVTensemble::Init()
{
	_position_controller->Init();
	_velocity_controller->Init();
	_force_controller->Init();


}

int NVTensemble::Run()
{

}