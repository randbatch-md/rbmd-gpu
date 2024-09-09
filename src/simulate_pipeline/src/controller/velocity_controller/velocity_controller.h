#pragma once
#include "model/device_data.h"
#include "model/structure_info_data.h"
#include <memory>
#include "data_manager.h"
#include "model/md_data.h"

class VelocityController
{
public:
	VelocityController()
	:_device_data(DataManager::getInstance().getDeviceData())//_device_data(std::make_shared<DeviceData>())		
	, _structure_info_data(DataManager::getInstance().getMDData()->_structure_info_data)//_structure_info_data(std::make_shared<StructureInfoData>())
	{};

	virtual ~VelocityController()=default;

	/**
	 * @brief Update Speed
	*/
	virtual void Update()=0;

	/**
	 * @brief Parameters and objects required for initializing the speed controller
	*/
	virtual void Init()=0;

protected:
	std::shared_ptr<StructureInfoData> _structure_info_data;
	std::shared_ptr<DeviceData> _device_data;

};