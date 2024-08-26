#pragma once
#include "basic_structure_data.h"

class FullStructureData : public BasicStructureData
{
public:
	bool checkStructure() const override
	{
		if (false == BasicStructureData::checkStructure())
		{
			return false;
		}

		//
		return true;
	}

	///bond
	rbmd::Id* _bond_types;
	rbmd::Id* _bond_first_id;
	rbmd::Id* _bond_second_id;

	///angle
	rbmd::Id* _angle_types;
	rbmd::Id* _angle_first_id;
	rbmd::Id* _angle_second_id;
	rbmd::Id* _angle_third_id;
};