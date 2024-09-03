#pragma once
#include <vector>
#include "common/types.h"
#include "common/object.h"

class ForceFieldData : public Object
{
public:
	virtual bool checkForceField() const = 0;

};
