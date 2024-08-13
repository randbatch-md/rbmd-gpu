#pragma once

#include "force_field_data.h"

class CVFFForceFieldData : public ForceFieldData
{
	bool checkForceField() const override
	{
		return true;
	}
};
