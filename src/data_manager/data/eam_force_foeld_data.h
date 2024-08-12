#pragma once

#include "force_field_data.h"

class EAMForceFieldData : public ForceFieldData
{
	bool checkForceField() const override
	{
		return true;
	}
};
