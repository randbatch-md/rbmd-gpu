#pragma once
#include "structure_data/structure_data.h"
#include "shake_controller.h"
#include "force_controller.h"
#include "force_field/force_field_data.h"
#include "model/device_data.h"
#include "model/structure_info_data.h"
#include <memory>

class VelocityController
{
public:
	VelocityController():_device_data(std::make_shared<DeviceData>()){};
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
	std::shared_ptr<StructureData> _structure_data;
	std::shared_ptr<ShakeController> _shake_controller;
	std::shared_ptr<ForceController> _force_controller;
	std::shared_ptr<ForceFieldData> _force_field_data;
	std::shared_ptr<StructureInfoData> _structure_info_data;

	std::shared_ptr<DeviceData> _device_data;
};