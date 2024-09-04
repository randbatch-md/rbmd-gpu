#pragma once
#include "force.h"
#include "common/types.h"

class LJForce: Force
{
public:
	LJForce();
	virtual ~LJForce() = default;

	void Init() override;
	void  Execute() override;

private:
	rbmd::Id _num_atoms;
};