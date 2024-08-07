#pragma once
#include "object.h"

class MemoryScheduler : public Object
{
public:
	virtual int AsyncMemoryH2D() = 0;
	virtual int AsyncMemoryD2H() = 0;
};