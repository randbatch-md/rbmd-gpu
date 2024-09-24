#pragma once
#include "../../common/types.h"
#include "force.h"
#include "model/box.h"
#include "neighbor_list/include/neighbor_list/neighbor_list.h"
#include "neighbor_list/include/neighbor_list_builder/full_neighbor_list_builder.h"
class LJForce : public Force
{
public:
	LJForce();
	virtual ~LJForce() = default;

	void Init() override;
	void  Execute() override;

private:
	rbmd::Id _num_atoms;
	std::shared_ptr<BaseNeighborListBuilder> _neighbor_list_builder;
	std::shared_ptr<NeighborList> list;
	Box box;

};