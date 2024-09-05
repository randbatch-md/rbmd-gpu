#pragma once
#include <vector>
#include "types.h"
#include "../object.h"

class ForceFieldData : public Object
{
public:
	virtual bool checkForceField() const = 0;

};
