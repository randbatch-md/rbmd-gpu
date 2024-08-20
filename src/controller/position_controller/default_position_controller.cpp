#include "default_position_controller.h"
#include "position_controller_op/update_position_op.h"

DefaultPositionController::DefaultPositionController() {};

void DefaultPositionController::Init()
{
	_dt = 0.001;           // �����ļ��ж�ȡ
	_num_atoms = _structure_info_data->_num_atoms;
}

void DefaultPositionController::Update() 
{
	// ����ʲôʱ����device��position �� velocity ʲôʱ���� host��
	bool available_shake = false; //�����ļ��ж�ȡ
	if (available_shake)
	{
		bool shake = false;// �����ļ��ж�ȡ
		if (shake)
		{
			// ��ǰlj����Ҫ����shake    
			//_device_data->_shake_vx = _device_data->_d_px;
			//_device_data->_shake_vy = _device_data->_d_py;
			//_device_data->_shake_vz = _device_data->_d_pz;
		}
		op::UpdatePositionOp<device::DEVICE_GPU> UpdatePositionOp;
		UpdatePositionOp<device::DEVICE_GPU>(_structure_info_data->_num_atoms,
			                                 _dt,
			                                 _structure_info_data->_range[0][0], // ����nvp�����Ǳ仯��
			                                 _structure_info_data->_range[1][0],
			                                 _structure_info_data->_range[2][0],
			                                 _structure_info_data->_range[0][1],
			                                 _structure_info_data->_range[1][1],
			                                 _structure_info_data->_range[2][1],
			                                 _device_data->_d_vx,
			                                 _device_data->_d_vy,
			                                 _device_data->_d_vz
			                                 _device_data->_d_px,
			                                 _device_data->_d_py,
			                                 _device_data->_d_pz);
	}
	else
	{
		op::UpdatePositionFlagOp<device::DEVICE_GPU> UpdatePositionFlagOp;
		UpdatePositionFlagOp<device::DEVICE_GPU>(_structure_info_data->_num_atoms,
			                                     _dt,
			                                     _range[0][0],
			                                     _range[1][0],
			                                     _range[2][0],
			                                     _range[0][1],
			                                     _range[1][1],
			                                     _range[2][1],
			                                     _device_data->_d_vx,
			                                     _device_data->_d_vy,
			                                     _device_data->_d_vz
			                                     _device_data->_d_px,
			                                     _device_data->_d_py,
			                                     _device_data->_d_pz,
			                                     _device_data->_d_flagX,// ������Ҫ����
			                                     _device_data->_d_flagY,
			                                     _device_data->_d_flagZ);
	}
}


void DefaultPositionController::SetCenterTargetPositions()
{
	std::string init_type = "inbuild";
	if (init_type==_init_type)
	{
		//Ϊ����rdf�����������׼����
	}
}
