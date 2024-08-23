#pragma once
#include "structure_data.h"
#include "shake_controller.h"
#include "device_data.h"
#include "structure_info_data.h"

class PositionController
{
public:
	PositionController() {};
	virtual ~PositionController()=default;

	virtual void Update() = 0;
	virtual void Init() = 0;

protected:
	std::shared_ptr<StructureData> _structure_data;
	std::shared_ptr<StructureInfoData> _structure_info_data;
	std::shared_ptr<ShakeController> _shake_controller;
	std::shared_ptr<DeviceData> _device_data;
};