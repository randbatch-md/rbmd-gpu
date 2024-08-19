#include "default_position_controller.h"
#include "position_controller_op/update_position_op.h"

DefaultPositionController::DefaultPositionController() {};
void DefaultPositionController::Update() 
{
	// 这里什么时候用device的position 和 velocity 什么时候用 host的
	bool available_shake = false; //配置文件中读取
	if (available_shake)
	{
		rbmd::Real min_x = 0;
		rbmd::Real min_y = 0;
		rbmd::Real min_z = 0;
		rbmd::Real max_x = 10;
		rbmd::Real max_y = 10;
		rbmd::Real max_z = 10;

		op::UpdatePositionFlagOp<device::DEVICE_GPU> UpdateVelocityOp;
		UpdatePositionFlagOp<device::DEVICE_GPU>(_structure_info_data->_num_atoms,
			                                     _dt,
			                                     _fmt2v,
												 _structure_info_data->_range[0][0],
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
			                                     _device_data->_d_pz,
			                                     flag_p);// 这里需要调用


		//auto&& position_flag = _para.GetFieldAsArrayHandle<Id3>(field::position_flag);
		//RunWorklet::UpdatePositionFlag(_dt, _velocity, _locator, _position, position_flag);
	}
	else
	{
		bool shake = false;// 配置文件中读取
		if (shake)
		{
			// 这些都是初始化的参数 看是否需要放到上一层的；
			auto& h_px = _structure_data->_h_px;
			auto& h_py = _structure_data->_h_py;
			auto& h_pz = _structure_data->_h_pz;

			// 可以直接传入 _force_controller ？    
			_shake_controller->_shake_vx = h_px;
			_shake_controller->_shake_vy = h_py;
			_shake_controller->_shake_vz = h_pz;
		}

		rbmd::Real min_x = 0;
		rbmd::Real min_y = 0;
		rbmd::Real min_z = 0;
		rbmd::Real max_x = 10;
		rbmd::Real max_y = 10;
		rbmd::Real max_z = 10;

		op::UpdatePositionOp<device::DEVICE_GPU> UpdateVelocityOp;
		UpdatePositionOp<device::DEVICE_GPU>(_structure_info_data->_num_atoms,
			                                 _dt,
			                                 _fmt2v,
											 min_x,
											 min_y,
											 min_z,
											 max_x,
											 max_y,
											 max_z,
			                                 _device_data->_d_vx,
			                                 _device_data->_d_vy,
			                                 _device_data->_d_vz
			                                 _device_data->_d_px,
			                                 _device_data->_d_py,
			                                 _device_data->_d_pz);

		//RunWorklet::UpdatePosition(_dt, _velocity, _locator, _position);
	}


}

void DefaultPositionController::Init()
{
	_dt = 0.001           // 配置文件中读取
}