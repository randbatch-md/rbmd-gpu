#pragma once
#include "common/device_types.h"
#include "common/types.h"

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
	struct UpdataVelocityRescaleOp
	{
		void operator()(const rbmd::Id& num_atoms,
			            const rbmd::Real& coeff_rescale,
			            rbmd::Real* vx,
			            rbmd::Real* vy,
			            rbmd::Real* vz);
	};

	template <typename DEVICE>
	struct UpdataVelocityNoseHooverOp
	{
		void operator()(const rbmd::Id& num_atoms,
			            const rbmd::Real& dt,
			            const rbmd::Real& fmt2v,
			            const rbmd::Real& nosehooverxi,
			            const rbmd::Real* mass,
						const rbmd::Real* fx,
						const rbmd::Real* fy,
						const rbmd::Real* fz,
			            rbmd::Real* vx,
			            rbmd::Real* vy,
			            rbmd::Real* vz);
	};

	template <typename DEVICE>
	struct UpdataVelocityBerendsenOp
	{
		void operator()(const rbmd::Id& num_atoms,
			            const rbmd::Real& coeff_Berendsen,
			            rbmd::Real* vx,
			            rbmd::Real* vy,
			            rbmd::Real* vz);
	};

}