#pragma once
#include "basic_structure_data.h"

class ChargeStructureData : public BasicStructureData
{
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