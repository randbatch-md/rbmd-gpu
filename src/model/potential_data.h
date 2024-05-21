#pragma once
#include <map>
#include "Types.h"

struct PotentialData
{
	bool hasMass()
	{
		return !_mass.empty();
	}


	std::map<rbmd::Id, rbmd::Real> _mass;
	std::map<rbmd::Id, rbmd::Real> _eps;
	std::map<rbmd::Id, rbmd::Real> _sigma;
};