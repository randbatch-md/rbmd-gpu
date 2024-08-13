#pragma once
#include "structure_data.h"

class BasicStructureData : public StructureData
{
public:
	bool checkStructure() const override
	{
		return true;
	}

	std::vector<rbmd::Real> _h_charge;

};
