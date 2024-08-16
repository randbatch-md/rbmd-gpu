#pragma once
#include "base/device_types.h"
#include "types.h"
#include "base/locator.h"
//#include "force_controller.h"
//#include <vector>
//#include <map>

namespace op
{
	template <typename DEVICE>
	struct UpdateVelocityOp
	{
		void operator()(const rbmd::Id& num_atoms,
			            const rbmd::Real& dt, 
			            const rbmd::Real& fmt2v,
			            const rbmd::Real* mass_map,
			            const rbmd::Real* fx,
		                const rbmd::Real* fy,
		                const rbmd::Real* fz,
			            rbmd::Real* vx,
			            rbmd::Real* vy,
			            rbmd::Real* vz);
	};
}