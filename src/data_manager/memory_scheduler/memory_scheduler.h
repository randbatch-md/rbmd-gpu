#pragma once
#include "object.h"

class MemoryScheduler : public Object
{
public:
	virtual int asyncMemoryH2D() = 0;
	virtual int asyncMemoryD2H() = 0;
};