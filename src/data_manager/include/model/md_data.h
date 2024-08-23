#pragma once
#include <memory>
#include "structure_data.h"
#include "structure_info_data.h"
#include "force_field_data.h"
#include "object.h"
#include "atoms_structure_data.h"
#include "charge_structure_data.h"
#include "full_structure_data.h"
#include "lj_force_foeld_data.h"
#include "cvff_force_foeld_data.h"
#include "eam_force_foeld_data.h"

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
		_force_field_data = std::make_shared<LJForceFieldData>();
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
			return false;
		}

		if (false == _force_field_data->checkForceField())
		{
			//log
			_console->error("check force field data failed!");
			return false;
		}

		return true;
	}

	std::shared_ptr<StructureData> _structure_data;
	std::shared_ptr<StructureInfoData> _structure_info_data;
	std::shared_ptr<ForceFieldData> _force_field_data;
};
