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
	void  ComputeChargeStructureFactorEwald(
		Box* box, 
		rbmd::Id _num_atoms,
		rbmd::Id Kmax, 
		rbmd::Real* value_Re_array,
		rbmd::Real* value_Im_array);

private:
	rbmd::Id _num_atoms;
	std::shared_ptr<BaseNeighborListBuilder> _neighbor_list_builder;
	std::shared_ptr<NeighborList> list;
	Box box;

	rbmd::Real _h_corr_value_x;
	rbmd::Real _h_corr_value_y;
	rbmd::Real _h_corr_value_z;

};