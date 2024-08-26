#pragma once
#include <hip/hip_runtime.h>
#include "object.h"

#include "include/data_manager.h"
#include "include/model/device_data.h"
#include "include/model/md_data.h"

class MemoryScheduler : public Object
{
public:
	/**
	 * @brief constructor
	*/
	MemoryScheduler();

	/**
	 * @brief async memeory host to device
	 * @return error code
	*/
	virtual bool asyncMemoryH2D() = 0;

	/**
	 * @brief async memory device to host
	 * @return error code
	*/
	virtual bool asyncMemoryD2H() = 0;

protected:
	DataManager& _data_manager;
	std::shared_ptr<DeviceData>& _device_data;
	const std::shared_ptr<MDData>& _md_data;
	const std::shared_ptr<StructureData>& _structure_data;
	const std::shared_ptr<StructureInfoData>& _structure_info_data;
	const std::shared_ptr<ForceFieldData>& _force_field_data;
};