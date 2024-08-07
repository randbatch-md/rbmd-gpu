#pragma once
#include "memory_scheduler.h"

class LJMemoryScheduler : public MemoryScheduler
{
public:
	int AsyncMemoryH2D() override;
	int AsyncMemoryD2H() override;
};
