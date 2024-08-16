#pragma once
#include "base/device_types.h"
#include "types.h"
#include "base/locator.h"
#include "force_controller.h"
#include <vector>
namespace op
{
	template <typename DEVICE>
	struct UpdateVelocityOp
	{
		void operator()(const rbmd::Real& dt, 
			            const rbmd::Real& fmt2v,
			            const std::vector<rbmd::Real>& fx,
		                const std::vector<rbmd::Real>& fy,
		                const std::vector<rbmd::Real>& fz,
			            const std::vector<rbmd::Real>& mass,
			            std::vector<rbmd::Real>& vx,
			            std::vector<rbmd::Real>& vy,
			            std::vector<rbmd::Real>& vz);
	};
}