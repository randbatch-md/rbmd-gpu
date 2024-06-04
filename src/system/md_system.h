#pragma once
#include "system.h"
#include "md_data.h"

class MDSystem : public System
{
public:
	MDSystem() = default;
	virtual ~MDSystem() = default;

public:
	int Evolve() override;
	auto& GetMDData() { return _md_data; }

private:
	int PreSolve();
	int Solve();
	int PostSolve();
private:
	MDData _md_data;
};