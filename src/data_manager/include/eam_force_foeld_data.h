#pragma once

#include "force_field_data.h"

class EAMForceFieldData : public ForceFieldData
{
public:

	bool checkForceField() const override
	{
		return true;
	}
};
