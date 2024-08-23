#pragma once
#include "device_types.h"
#include "types.h"

namespace op
{
	template <typename DEVICE>
	struct ComputeTemperatureOp
	{
		void operator()(const rbmd::Id& num_atoms,
			            const rbmd::Real& mvv2e,
			            const rbmd::Real* mass,
			            const rbmd::Real* vx,
		                const rbmd::Real* vy,
		                const rbmd::Real* vz,
			            rbmd::Real& temp_sum);
	};

	template <typename DEVICE>
	struct UpdataVelocityOp
	{
		void operator()(const rbmd::Id& num_atoms,
			            const rbmd::Real& coeff_rescale,
			            rbmd::Real* vx,
			            rbmd::Real* vy,
			            rbmd::Real* vz);
	};
}