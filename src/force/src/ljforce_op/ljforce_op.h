#pragma once
#include "common/types.h"


namespace op
{

	template <typename DEVICE>
	struct LJforceOp
	{
		void operator()(Box* box,
			            const rbmd::Id& N,
			            const rbmd::Id* atoms_type,
			            const rbmd::Id* molecular_type,
			            const rbmd::Real* sigma,
			            const rbmd::Real* eps,
			            const rbmd::Id* start_id,
		                const rbmd::Id* end_id,
		                const rbmd::Id* id_verletlist,
			            const rbmd::Real* px,
			            const rbmd::Real* py,
			            const rbmd::Real* pz,
			            rbmd::Real* force_x,
			            rbmd::Real* force_y,
			            rbmd::Real* force_z,
			            rbmd::Real* evdwl);
	};

    template <typename T>
    void DeviceReduceSum(const T* d_data, T* d_result, int num_elements);

}