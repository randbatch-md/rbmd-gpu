#pragma once
#include "../../common/device_types.h"
#include "../../common/types.h"
#include "../../data_manager/include/model/box.h"

namespace op
{

	template <typename DEVICE>
	struct LJForceOp
	{
		void operator()(Box* box,
			            const rbmd::Real cut_off,
			            const rbmd::Id num_atoms,
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
			            rbmd::Real* evdwl,
			            rbmd::Real* total_evdwl);
	};

	template <typename DEVICE>
	struct LJRBLForceOp
	{
		void operator()(Box* box,
			const rbmd::Real rs,
			const rbmd::Real rc,
			const rbmd::Id pice_num,
			const rbmd::Id num_atoms,
			const rbmd::Id* atoms_type,
			const rbmd::Id* molecular_type,
			const rbmd::Real* sigma,
			const rbmd::Real* eps,
			const rbmd::Id* start_id,
			const rbmd::Id* end_id,
			const rbmd::Id* id_verletlist,
			const rbmd::Id* id_verletlist_rcs,
			const rbmd::Real* px,
			const rbmd::Real* py,
			const rbmd::Real* pz,
			rbmd::Real* corr_force_x,
			rbmd::Real* corr_force_y,
			rbmd::Real* corr_force_z,
			rbmd::Real* corr_value_x,
			rbmd::Real* corr_value_y,
			rbmd::Real* corr_value_z);
	};

	template <>
	struct LJForceOp<device::DEVICE_GPU>
	{
		void operator()(Box* box,
			const rbmd::Real cut_off,
			const rbmd::Id num_atoms,
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
			rbmd::Real* evdwl,
			rbmd::Real* total_evdwl);
	};

	template <>
	struct LJRBLForceOp<device::DEVICE_GPU>
	{
		void operator()(Box* box,
			const rbmd::Real rs,
			const rbmd::Real rc,
			const rbmd::Id pice_num,
			const rbmd::Id num_atoms,
			const rbmd::Id* atoms_type,
			const rbmd::Id* molecular_type,
			const rbmd::Real* sigma,
			const rbmd::Real* eps,
			const rbmd::Id* start_id,
			const rbmd::Id* end_id,
			const rbmd::Id* id_verletlist,
			const rbmd::Id* id_verletlist_rcs,
			const rbmd::Real* px,
			const rbmd::Real* py,
			const rbmd::Real* pz,
			rbmd::Real* corr_force_x,
			rbmd::Real* corr_force_y,
			rbmd::Real* corr_force_z,
			rbmd::Real* corr_value_x,
			rbmd::Real* corr_value_y,
			rbmd::Real* corr_value_z);
	};
    //template <typename T>
    //void DeviceReduceSum(const T* d_data, T* d_result, int num_elements);

}// namespace op