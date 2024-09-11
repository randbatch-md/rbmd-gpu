#pragma once
#include "force.h"
#include "../../common/types.h"
#include "../neighbor_list/include/neighbor_list.h"
#include "../neighbor_list/include/neighbor_list_builder/full_neighbor_list_builder.h"
#include "model/box.h"

class LJForce : public Force
{
public:
	LJForce();
	virtual ~LJForce() = default;

	void Init() override;
	void  Execute() override;

private:
	rbmd::Id _num_atoms;
	std::shared_ptr<FullNeighborListBuilder> full_list_builder;
	std::shared_ptr<NeighborList> list;
	Box box;

};