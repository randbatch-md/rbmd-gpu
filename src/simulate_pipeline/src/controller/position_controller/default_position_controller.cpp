#include "default_position_controller.h"
#include "update_position_op.h"
#include <thrust/device_ptr.h>
#define POSITION_OUTPUT
DefaultPositionController::DefaultPositionController() {};

void DefaultPositionController::Init()
{
	//auto atom_style = _config_data->Get<std::string>("atom_style", "init_configuration", "read_data");

	_dt = 0.001;           // 配置文件中读取
	_num_atoms = *(_structure_info_data->_num_atoms);
}

void DefaultPositionController::Update() 
{
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
		op::UpdatePositionOp<device::DEVICE_GPU> update_position_op;
		update_position_op(_num_atoms,
			               _dt,
			               (*_structure_info_data->_range)[0][0], // 这在nvp里面是变化的
			               (*_structure_info_data->_range)[1][0],
			               (*_structure_info_data->_range)[2][0],
			               (*_structure_info_data->_range)[0][1],
			               (*_structure_info_data->_range)[1][1],
			               (*_structure_info_data->_range)[2][1],
			               thrust::raw_pointer_cast(_device_data->_d_vx.data()),
			               thrust::raw_pointer_cast(_device_data->_d_vy.data()),
			               thrust::raw_pointer_cast(_device_data->_d_vz.data()),
			               thrust::raw_pointer_cast(_device_data->_d_px.data()),
			               thrust::raw_pointer_cast(_device_data->_d_py.data()),
			               thrust::raw_pointer_cast(_device_data->_d_pz.data()));
	}
	else
	{
		op::UpdatePositionFlagOp<device::DEVICE_GPU> update_position_op;
		update_position_op(_num_atoms,
			               _dt,
			               (*_structure_info_data->_range)[0][0],
			               (*_structure_info_data->_range)[1][0],
			               (*_structure_info_data->_range)[2][0],
			               (*_structure_info_data->_range)[0][1],
			               (*_structure_info_data->_range)[1][1],
			               (*_structure_info_data->_range)[2][1],
			               thrust::raw_pointer_cast(_device_data->_d_vx.data()),
			               thrust::raw_pointer_cast(_device_data->_d_vy.data()),
			               thrust::raw_pointer_cast(_device_data->_d_vz.data()),
			               thrust::raw_pointer_cast(_device_data->_d_px.data()),
			               thrust::raw_pointer_cast(_device_data->_d_py.data()),
			               thrust::raw_pointer_cast(_device_data->_d_pz.data()),
			               thrust::raw_pointer_cast(_device_data->_d_flagX.data()),// 这里需要调用
			               thrust::raw_pointer_cast(_device_data->_d_flagY.data()),
			               thrust::raw_pointer_cast(_device_data->_d_flagZ.data()));
	}

#ifdef POSITION_OUTPUT
	// 分配主机内存
	rbmd::Real* h_vx = new rbmd::Real[_device_data->_d_vx.size()];
	rbmd::Real* h_vy = new rbmd::Real[_device_data->_d_vy.size()];
	rbmd::Real* h_vz = new rbmd::Real[_device_data->_d_vz.size()];

	rbmd::Real* h_px = new rbmd::Real[_device_data->_d_px.size()];
	rbmd::Real* h_py = new rbmd::Real[_device_data->_d_py.size()];
	rbmd::Real* h_pz = new rbmd::Real[_device_data->_d_pz.size()];

	// 将设备数据拷贝到主机
	CHECK_RUNTIME(MEMCPY(h_vx, thrust::raw_pointer_cast(_device_data->_d_vx.data()), _device_data->_d_vx.size() * sizeof(rbmd::Real), D2H));
	CHECK_RUNTIME(MEMCPY(h_vy, thrust::raw_pointer_cast(_device_data->_d_vy.data()), _device_data->_d_vy.size() * sizeof(rbmd::Real), D2H));
	CHECK_RUNTIME(MEMCPY(h_vz, thrust::raw_pointer_cast(_device_data->_d_vz.data()), _device_data->_d_vz.size() * sizeof(rbmd::Real), D2H));

	CHECK_RUNTIME(MEMCPY(h_px, thrust::raw_pointer_cast(_device_data->_d_px.data()), _device_data->_d_px.size() * sizeof(rbmd::Real), D2H));
	CHECK_RUNTIME(MEMCPY(h_py, thrust::raw_pointer_cast(_device_data->_d_py.data()), _device_data->_d_py.size() * sizeof(rbmd::Real), D2H));
	CHECK_RUNTIME(MEMCPY(h_pz, thrust::raw_pointer_cast(_device_data->_d_pz.data()), _device_data->_d_pz.size() * sizeof(rbmd::Real), D2H));

	std::ofstream output_file_p_velocity;
	output_file_p_velocity.open("output_file_p_velocity_" + std::to_string(test_current_step) + ".txt");
	for (size_t i = 0; i < _num_atoms; i++)
	{
		output_file_p_velocity << i << "," << h_vx[i] << "," << h_vy[i] << "," << h_vz[i] << std::endl;

	}

	std::ofstream output_file_p_position;
	output_file_p_position.open("output_file_p_position_" + std::to_string(test_current_step) + ".txt");
	for (size_t i = 0; i < _num_atoms; i++)
	{
		output_file_p_position << i << "," << h_px[i] << "," << h_py[i] << "," << h_pz[i] << std::endl;

	}

#endif // V_OUTPUT

}


void DefaultPositionController::SetCenterTargetPositions()
{
	std::string init_type = "inbuild";
	if (init_type==_init_type)
	{
		//为计算rdf做填充数据做准备；
	}
}
