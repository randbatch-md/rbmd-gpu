#pragma once
#include <vector>
#include "types.h"
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
	rbmd::Real* _h_mass;
	rbmd::Real* _h_eps;
	rbmd::Real* _h_sigma;
};
