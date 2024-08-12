#pragma once
#include <map>
#include "../model/types.h"

struct ForceFieldData
{
	bool hasMass()
	{
		return !_h_mass.empty();
	}

	std::map<rbmd::Id, rbmd::Real> _h_mass;
	std::map<rbmd::Id, rbmd::Real> _h_eps;
	std::map<rbmd::Id, rbmd::Real> _h_sigma;
};
