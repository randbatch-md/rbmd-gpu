#pragma once

#include "force_field_data.h"

class CVFFForceFieldData : public ForceFieldData
{
public:
	bool checkForceField() const override
	{
		return true;
	}

	///bond
	rbmd::Real* _bond_coeffs_;
	rbmd::Real* _bond_coeffs_equilibrium;

	///angle
	rbmd::Real* _angle_coeffs_k;
	rbmd::Real* _angle_coeffs_equilibrium;
};
