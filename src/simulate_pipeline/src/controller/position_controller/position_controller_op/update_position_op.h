#pragma once
#include "device_types.h"
#include "types.h"

namespace op
{
	template <typename DEVICE>
	struct UpdatePositionFlagOp
	{
		void operator()(const rbmd::Id num_atoms,
			            const rbmd::Real dt,
			            const rbmd::Real min_x,
			            const rbmd::Real min_y,
			            const rbmd::Real min_z,
			            const rbmd::Real max_x,
			            const rbmd::Real max_y,
			            const rbmd::Real max_z,
			            const rbmd::Real* vx,
			            const rbmd::Real* vy,
			            const rbmd::Real* vz,
		                rbmd::Real* px,
		                rbmd::Real* py,
		                rbmd::Real* pz,
			            rbmd::Id* flag_px,
			            rbmd::Id* flag_py,
			            rbmd::Id* flag_pz);
	};

	template <typename DEVICE>
	struct UpdatePositionOp
	{
		void operator()(const rbmd::Id num_atoms,
			            const rbmd::Real dt,
			            const rbmd::Real min_x,
			            const rbmd::Real min_y,
			            const rbmd::Real min_z,
			            const rbmd::Real max_x,
			            const rbmd::Real max_y,
			            const rbmd::Real max_z,
			            const rbmd::Real* vx,
			            const rbmd::Real* vy,
			            const rbmd::Real* vz,
			            rbmd::Real* px,
			            rbmd::Real* py,
			            rbmd::Real* pz);
	};

	template <>
	struct UpdatePositionFlagOp<device::DEVICE_GPU>
	{
		void operator()(const rbmd::Id num_atoms,
						const rbmd::Real dt,
						const rbmd::Real min_x,
						const rbmd::Real min_y,
						const rbmd::Real min_z,
						const rbmd::Real max_x,
						const rbmd::Real max_y,
						const rbmd::Real max_z,
						const rbmd::Real* vx,
						const rbmd::Real* vy,
						const rbmd::Real* vz,
						rbmd::Real* px,
						rbmd::Real* py,
						rbmd::Real* pz,
						rbmd::Id* flag_px,
						rbmd::Id* flag_py,
						rbmd::Id* flag_pz);
	};

	template <>
	struct UpdatePositionOp<device::DEVICE_GPU>
	{
		void operator()(const rbmd::Id num_atoms,
						const rbmd::Real dt,
						const rbmd::Real min_x,
						const rbmd::Real min_y,
						const rbmd::Real min_z,
						const rbmd::Real max_x,
						const rbmd::Real max_y,
						const rbmd::Real max_z,
						const rbmd::Real* vx,
						const rbmd::Real* vy,
						const rbmd::Real* vz,
						rbmd::Real* px,
						rbmd::Real* py,
						rbmd::Real* pz);
	};
}