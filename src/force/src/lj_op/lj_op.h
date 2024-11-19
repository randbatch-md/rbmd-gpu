#pragma once
#include "../../common/device_types.h"
#include "../../common/types.h"
#include "../../data_manager/include/model/box.h"
#include "../common/erf_table.h"

namespace op
{

        template <typename DEVICE>
        struct LJForceOp
        {
          void operator()( Box box,
                          const rbmd::Real cut_off,
                          const rbmd::Id num_atoms,
                          const rbmd::Id* atoms_type,
                          const rbmd::Real* sigma,
                          const rbmd::Real* eps,
                          const rbmd::Id* start_id,
                          const rbmd::Id* end_id,
                          const rbmd::Id* id_verletlist,
                          const rbmd::Real* px,
                          const rbmd::Real* py,
                          const rbmd::Real* pz,
                          rbmd::Real* fx,
                          rbmd::Real* fy,
                          rbmd::Real* fz,
                          rbmd::Real* flat_virial,
                          rbmd::Real* total_evdwl);
        };

	template <typename DEVICE>
	struct LJRBLForceOp
	{
		void operator()( Box box,
			const rbmd::Real rs,
			const rbmd::Real rc,
			const rbmd::Id num_atoms,
			const rbmd::Id neighbor_sample_num,
			const rbmd::Id pice_num,
			const rbmd::Id* atoms_type,
			const rbmd::Real* sigma,
			const rbmd::Real* eps,
			const rbmd::Id* start_id,
			const rbmd::Id* end_id,
			const rbmd::Id* id_verletlist,
			const rbmd::Id* id_random_neighbor,
			const rbmd::Id* random_neighbor_num,
			const rbmd::Real* px,
			const rbmd::Real* py,
			const rbmd::Real* pz,
			rbmd::Real* fx,
                        rbmd::Real* fy,
                        rbmd::Real* fz);
	};

	template <typename DEVICE>
	struct FixRBLForceOp
	{
		void operator()(
			const rbmd::Id num_atoms,
			const rbmd::Real corr_value_x,
			const rbmd::Real corr_value_y,
			const rbmd::Real corr_value_z,
			rbmd::Real* fx,
                        rbmd::Real* fy,
                        rbmd::Real* fz);
	};

	template <typename DEVICE>
	struct LJEnergyOp
	{
		void operator()( Box box,
			const rbmd::Real cut_off,
			const rbmd::Id num_atoms,
			const rbmd::Id* atoms_type,
			const rbmd::Real* sigma,
			const rbmd::Real* eps,
			const rbmd::Id* start_id,
			const rbmd::Id* end_id,
			const rbmd::Id* id_verletlist,
			const rbmd::Real* px,
			const rbmd::Real* py,
			const rbmd::Real* pz,
			rbmd::Real* flat_virial,
			rbmd::Real* total_evdwl);
	};

        // // // // // // // // // // // // // // // // // // //
        template <>
        struct LJForceOp<device::DEVICE_GPU>
        {
          void operator()(const Box  box,
                          const rbmd::Real cut_off,
                          const rbmd::Id num_atoms,
                          const rbmd::Id* atoms_type,
                          const rbmd::Real* sigma,
                          const rbmd::Real* eps,
                          const rbmd::Id* start_id,
                          const rbmd::Id* end_id,
                          const rbmd::Id* id_verletlist,
                          const rbmd::Real* px,
                          const rbmd::Real* py,
                          const rbmd::Real* pz,
                          rbmd::Real* fx,
                          rbmd::Real* fy,
                          rbmd::Real* fz,
                          rbmd::Real* flat_virial,
                          rbmd::Real* total_evdwl);
        };

	template <>
	struct LJRBLForceOp<device::DEVICE_GPU>
	{
		void operator()(const  Box box,
			const rbmd::Real rs,
			const rbmd::Real rc,
			const rbmd::Id num_atoms,
			const rbmd::Id neighbor_sample_num,
			const rbmd::Id pice_num,
			const rbmd::Id* atoms_type,
			const rbmd::Real* sigma,
			const rbmd::Real* eps,
			const rbmd::Id* start_id,
			const rbmd::Id* end_id,
			const rbmd::Id* id_verletlist,
			const rbmd::Id* id_random_neighbor,
			const rbmd::Id* random_neighbor_num,
			const rbmd::Real* px,
			const rbmd::Real* py,
			const rbmd::Real* pz,
			rbmd::Real* fx,
                        rbmd::Real* fy,
                        rbmd::Real* fz);
	};

	template <>
	struct FixRBLForceOp<device::DEVICE_GPU>
	{
		void operator()(
			const rbmd::Id num_atoms,
			const rbmd::Real corr_value_x,
			const rbmd::Real corr_value_y,
			const rbmd::Real corr_value_z,
			rbmd::Real* fx,
                        rbmd::Real* fy,
                        rbmd::Real* fz);
	};

	template <>
	struct LJEnergyOp<device::DEVICE_GPU>
	{
		void operator()( Box box,
			const rbmd::Real cut_off,
			const rbmd::Id num_atoms,
			const rbmd::Id* atoms_type,
			const rbmd::Real* sigma,
			const rbmd::Real* eps,
			const rbmd::Id* start_id,
			const rbmd::Id* end_id,
			const rbmd::Id* id_verletlist,
			const rbmd::Real* px,
			const rbmd::Real* py,
			const rbmd::Real* pz,
			rbmd::Real* flat_virial,
			rbmd::Real* total_evdwl);
	};

}// namespace op