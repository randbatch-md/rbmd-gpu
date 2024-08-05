#pragma once
#include <vector>
#include "types.h"

struct PotentialData
{
	bool hasMass()
	{
		return !_mass.empty();
	}


	std::vector<rbmd::Real> _mass;
	std::vector<rbmd::Real> _eps;
	std::vector<rbmd::Real> _sigma;
};