#include "./include/scheduler/eam_memory_scheduler.h"

bool EAMMemoryScheduler::asyncMemoryH2D()
{
	if (false == MemoryScheduler::asyncMemoryH2D())
	{
		//log
		return false;
	}
	auto fd = std::dynamic_pointer_cast<EAMForceFieldData>(_force_field_data);

	///copy force field
	_device_data->_d_frho.resize(fd->_num_frho);
	_device_data->_d_rhor.resize(fd->_num_rhor);
	_device_data->_d_z2r.resize(fd->_num_z2r);

	///frho
	thrust::copy(fd->_h_frho, fd->_h_frho + fd->_num_frho, _device_data->_d_frho.begin());

	///rhor
	thrust::copy(fd->_h_rhor, fd->_h_rhor + fd->_num_rhor, _device_data->_d_rhor.begin());

	///z2r
	thrust::copy(fd->_h_z2r, fd->_h_z2r + fd->_num_z2r, _device_data->_d_z2r.begin());

	return true;
}

bool EAMMemoryScheduler::asyncMemoryD2H()
{
	return true;
}
