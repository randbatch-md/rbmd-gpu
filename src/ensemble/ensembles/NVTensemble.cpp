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

	// �ٶȺ�λ�õĳ�ʼ�� ������ _positiong  velocity�ĳ�ʼ�����Լ��໥֮��ĳ�ʼ��������charge��mass��atoms_id��_molecule_id��_Volume �� ���������Ǹ�����Ҫ�Ĳ���
	_position_controller.Init(); 
	_velocity_controller.Init();

	// ���ĳ�ʼ�� ������ͬ�������ļ���ȡ(��EAM�����ļ�)����Զ�����ĳ�ʼ�������ͼ��㣨��Զ�̣�_psamplekey��
	_force_controller.Init(); 


}

int NVTensemble::Run()
{

}