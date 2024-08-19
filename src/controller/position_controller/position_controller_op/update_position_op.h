#pragma once
#include "base/device_types.h"
#include "types.h"
#include "base/locator.h"

namespace op
{
	template <typename DEVICE>
	struct UpdatePositionFlagOp
	{
		void operator()(const rbmd::Id& num_atoms,
			            const rbmd::Real& dt,
			            const rbmd::Real* min_x,
			            const rbmd::Real* min_y,
			            const rbmd::Real* min_z,
			            const rbmd::Real* max_x,
			            const rbmd::Real* max_y,
			            const rbmd::Real* max_z,
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
		void operator()(const rbmd::Id& num_atoms,
			            const rbmd::Real& dt,
			            const rbmd::Real* min_x, // 没有在设备上分配内存的使用的是引用还是指针
			            const rbmd::Real* min_y,
			            const rbmd::Real* min_z,
			            const rbmd::Real* max_x,
			            const rbmd::Real* max_y,
			            const rbmd::Real* max_z,
			            const rbmd::Real* vx,
			            const rbmd::Real* vy,
			            const rbmd::Real* vz,
			            rbmd::Real* px,
			            rbmd::Real* py,
			            rbmd::Real* pz);
	};
}