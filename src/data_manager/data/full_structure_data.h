#pragma once
#include "basic_structure_data.h"

class FullStructureData : public BasicStructureData
{
public:
	bool checkData() const override
	{
		if (false == BasicStructureData::checkData())
		{
			return false
		}

		//
		return true;
	}
};