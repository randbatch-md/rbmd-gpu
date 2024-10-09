#pragma once

#include "force_field_data.h"
#include "common/types.h"

class LJForceFieldData : public ForceFieldData
{
public:

	bool checkForceField() const override
	{

		return true;
	}

public:
	///mass
	rbmd::Real* _h_mass;

	///eps
	rbmd::Real* _h_eps;

	///sigma
	rbmd::Real* _h_sigma;

	//charge
	//rbmd::Real* _h_charge;
};
