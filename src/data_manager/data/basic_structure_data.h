#pragma once
#include "structure_data.h"

class BasicStructureData : public StructureData
{
public:
	bool checkData() const override
	{

	}

	std::vector<rbmd::Real> _h_charge;

};
