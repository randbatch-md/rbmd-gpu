#pragma once
#include "structure_data/structure_data.h"
#include "shake_controller.h"
#include "model/device_data.h"
#include "model/structure_info_data.h"

class PositionController
{
public:
	PositionController() {};
	virtual ~PositionController()=default;

	/**
	 * @brief Update Position
	*/
	virtual void Update() = 0;

	/**
	 * @brief Parameters and objects required for initializing the Position controller
	*/
	virtual void Init() = 0;

protected:
	std::shared_ptr<StructureData> _structure_data;
	std::shared_ptr<StructureInfoData> _structure_info_data;
	std::shared_ptr<ShakeController> _shake_controller;
	std::shared_ptr<DeviceData> _device_data;
};