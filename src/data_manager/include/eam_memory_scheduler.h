#pragma once
#include "memory_scheduler.h"

class EAMMemoryScheduler : public MemoryScheduler
{
public:
	int asyncMemoryH2D() override;
	int asyncMemoryD2H() override;
};
