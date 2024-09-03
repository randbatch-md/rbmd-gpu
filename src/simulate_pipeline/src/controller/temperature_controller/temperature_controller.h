#pragma once
#include "../data_manager/include/model/device_data.h"
#include "../data_manager/include/model/structure_info_data.h"
#include <memory>

class TemperatureController
{
public:
	TemperatureController() :_device_data(std::make_shared<DeviceData>()) {};
	virtual ~TemperatureController()=default;

	/**
	 * @brief Update Temperature
	*/
	virtual void Update()=0;

	/**
	 * @brief Parameters and objects required for initializing the temperature controller
	*/
	virtual void Init()=0;

protected:
	std::shared_ptr<StructureInfoData> _structure_info_data;
	std::shared_ptr<DeviceData> _device_data;
};