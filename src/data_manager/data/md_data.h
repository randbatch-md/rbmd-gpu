#pragma once
#include <memory>
#include "structure_data.h"
#include "structure_info_data.h"
#include "force_field_data.h"
#include "object.h"

class MDData : public Object
{
public:
	MDData()
	{
		_structure_data = std::make_shared<LJStructureData>();
		_structure_info_data = std::make_shared<StructureInfoData>();
		_force_field_data = std::make_shared<ForceFieldData>();
	}

	bool checkMDData()
	{
		if (0 == _structure_data->checkData())
		{
			//log
			_console->error("check structure data failed!");
		}

		if (0 == _force_field_data->checkForceField())
		{
			//log
			_console->error("check force field data failed!");
		}
	}

	std::shared_ptr<StructureData> _structure_data;
	std::shared_ptr<StructureInfoData> _structure_info_data;
	std::shared_ptr<ForceFieldData> _force_field_data;
};
