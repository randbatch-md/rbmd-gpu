#include "update_velocity_op.h"
#include "device_types.h"
#include<iostream>

namespace op
{		
	struct UpdateVelocityOp<device::DEVICE_CPU>
	{
		void operator()(const rbmd::Real& dt,
			const rbmd::Real& fmt2v,
			const std::vector<rbmd::Real>& fx,
			const std::vector<rbmd::Real>& fy,
			const std::vector<rbmd::Real>& fz,
			const std::vector<rbmd::Real>& mass,
			std::vector<rbmd::Real>& vx,
			std::vector<rbmd::Real>& vy,
			std::vector<rbmd::Real>& vz)
		{
			auto num_atoms = fx.size();
			for (size_t i = 0; i < num_atoms; i++)
			{
				vx += 0.5 * fx / mass * dt * fmt2v;
				vy += 0.5 * fy / mass * dt * fmt2v;
				vz += 0.5 * fz / mass * dt * fmt2v;
			}
		}
	};
	
	template struct UpdateVelocityOp<device::DEVICE_CPU>;
	template struct UpdateVelocityOp<device::DEVICE_CPU>;
}