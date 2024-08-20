#include "default_position_controller.h"
#include "position_controller_op/update_position_op.h"

DefaultPositionController::DefaultPositionController() {};

void DefaultPositionController::Init()
{
	_dt = 0.001;           // 配置文件中读取
	_num_atoms = _structure_info_data->_num_atoms;
}

void DefaultPositionController::Update() 
{
	// 这里什么时候用device的position 和 velocity 什么时候用 host的
	bool available_shake = false; //配置文件中读取
	if (available_shake)
	{
		bool shake = false;// 配置文件中读取
		if (shake)
		{
			// 当前lj不需要带入shake    
			//_device_data->_shake_vx = _device_data->_d_px;
			//_device_data->_shake_vy = _device_data->_d_py;
			//_device_data->_shake_vz = _device_data->_d_pz;
		}
		op::UpdatePositionOp<device::DEVICE_GPU> UpdatePositionOp;
		UpdatePositionOp<device::DEVICE_GPU>(_structure_info_data->_num_atoms,
			                                 _dt,
			                                 _structure_info_data->_range[0][0], // 这在nvp里面是变化的
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
			                                     _device_data->_d_flagX,// 这里需要调用
			                                     _device_data->_d_flagY,
			                                     _device_data->_d_flagZ);
	}
}


void DefaultPositionController::SetCenterTargetPositions()
{
	std::string init_type = "inbuild";
	if (init_type==_init_type)
	{
		//为计算rdf做填充数据做准备；
	}
}
