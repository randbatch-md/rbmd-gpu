#pragma once
#include "structure_data/structure_data.h"
#include "model/device_data.h"
#include "model/structure_info_data.h"
#include "force/include/box1.h"

class Force
{
public:
	Force() {};
	virtual ~Force()=default;

	//virtual void Update()=0;
	virtual void Init() {};
	virtual void Execute() = 0;

protected:
	std::shared_ptr<StructureData> _structure_data;
	std::shared_ptr<StructureInfoData> _structure_info_data;
	std::shared_ptr<DeviceData> _device_data;
	BOX box;

};