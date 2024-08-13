#pragma once
#include <map>
#include "../model/types.h"
#include "object.h"

class ForceFieldData : public Object
{
public:
	virtual bool checkForceField() const = 0;

private:
	bool hasMass()
	{
		return !_h_mass.empty();
	}

public:
	std::map<rbmd::Id, rbmd::Real> _h_mass;
	std::map<rbmd::Id, rbmd::Real> _h_eps;
	std::map<rbmd::Id, rbmd::Real> _h_sigma;
};
