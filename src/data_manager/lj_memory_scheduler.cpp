#include "include/lj_memory_scheduler.h"
#include "include/data_manager.h"
#include "include/device_data.h"
#include "include/md_data.h"

bool LJMemoryScheduler::asyncMemoryH2D()
{
	auto& data_manager = DataManager::getInstance();
	auto& device_data = data_manager.getDeviceData();
	auto& md_data = data_manager.getMDData();
	auto& num_atoms = md_data->_structure_info_data->_num_atoms;
	auto& structure_data = md_data->_structure_data;
	auto& h_px = structure_data->_h_px;

	device_data->_d_px.resize(num_atoms);
	thrust::copy(h_px, h_px + num_atoms, device_data->_d_px);

	return true;
}

bool LJMemoryScheduler::asyncMemoryD2H()
{
	return true;
}
