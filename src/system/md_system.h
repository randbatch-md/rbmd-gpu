#pragma once
#include "system.h"
#include "md_data.h"

class MDSystem : public System
{
public:
	MDSystem() = default;
	virtual ~MDSystem() = default;

	int Evolve() override;

public:
	auto& GetMDData() { return _md_data; }
private:
	MDData _md_data;
};