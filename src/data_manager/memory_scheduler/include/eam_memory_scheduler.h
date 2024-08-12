#pragma once
#include "memory_scheduler.h"

class EAMMemoryScheduler : public MemoryScheduler
{
public:
	int AsyncMemoryH2D() override;
	int AsyncMemoryD2H() override;
};
