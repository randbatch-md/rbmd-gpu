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
};