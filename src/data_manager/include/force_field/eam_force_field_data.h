#pragma once

#include "force_field_data.h"
#include "common/types.h"

class EAMForceFieldData : public ForceFieldData
{
public:

	bool checkForceField() const override
	{
		return true;
	}

	///F(ρ) on host
	rbmd::Id _num_frho;
	rbmd::Real* _h_frho;

	///ρ(r) on host
	rbmd::Id _num_rhor;
	rbmd::Real* _h_rhor;

	///ϕ(r) on host
	rbmd::Id _num_z2r;
	rbmd::Real* _h_z2r;
};
