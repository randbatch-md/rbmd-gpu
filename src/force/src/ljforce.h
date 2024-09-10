#pragma once
#include "force.h"
#include "../../common/types.h"
#include "neighbor_list/include/neighbor_list.h"



class LJForce: Force
{
public:
	LJForce();
	virtual ~LJForce() = default;

	void Init() override;
	void  Execute() override;

private:
	rbmd::Id _num_atoms;
	NeighborList list;
	Box box;

};