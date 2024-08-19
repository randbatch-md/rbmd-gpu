#pragma once

#include "force_field_data.h"

class LJForceFieldData : public ForceFieldData
{
public:

	bool checkForceField() const override
	{

		return true;
	}
};
