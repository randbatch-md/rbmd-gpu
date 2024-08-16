#pragma once
#include "memory_scheduler.h"

class CVFFMemoryScheduler : public MemoryScheduler
{
public:
	bool asyncMemoryH2D() override;
	bool asyncMemoryD2H() override;
};
