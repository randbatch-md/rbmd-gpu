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
	if (_para.GetParameter<std::string>(PARA_FILE_TYPE) == "EAM")
	{
		_potential_file.open(_para.GetParameter<std::string>(PARA_POTENTIAL_FILE));
		ReadPotentialFile(_potential_file);
		InitStyle();
	}

	ExecutionMD::Init();

	InitialCondition();

	ComputeForce();

	// 速度和位置的初始化 包括了 _positiong  velocity的初始化，以及相互之间的初始化；还有charge、mass、atoms_id、_molecule_id、_Volume 等 计算他们是各自需要的参数
	_position_controller.Init(); 
	_velocity_controller.Init();

	// 力的初始化 包括不同力场的文件读取(如EAM势能文件)、近远程力的初始参数填充和计算（如远程：_psamplekey）
	_force_controller.Init(); 


}

int NVTensemble::Run()
{

}