#pragma once
#include <memory>
#include "structure_data.h"
#include "structure_info_data.h"
#include "force_field_data.h"
#include "object.h"

class AtomsStructureData;
class ChargeStructureData;
class FullStructureData;

class MDData : public Object
{
public:
	/**
	 * @brief constructor
	*/
	MDData()
	{
		_structure_data = std::make_shared<AtomsStructureData>();
		_structure_info_data = std::make_shared<StructureInfoData>();
		_force_field_data = std::make_shared<ForceFieldData>();
	}

	/**
	 * @brief check data
	 * @return true or false
	*/
	bool checkMDData()
	{
		if (false == _structure_data->checkStructure())
		{
			//log
			_console->error("check structure data failed!");
		}

		if (false == _force_field_data->checkForceField())
		{
			//log
			_console->error("check force field data failed!");
		}
	}

	std::shared_ptr<StructureData> _structure_data;
	std::shared_ptr<StructureInfoData> _structure_info_data;
	std::shared_ptr<ForceFieldData> _force_field_data;
};
