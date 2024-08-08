#pragma once
#include <map>
#include "../model/types.h"

struct ForceFieldData
{
	std::map<rbmd::Id, rbmd::Real> _h_mass;
};
