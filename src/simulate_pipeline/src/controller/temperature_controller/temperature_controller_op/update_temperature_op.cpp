//#include "update_temperature_op.h"
//#include "device_types.h"
//#include<iostream>
//
//namespace op
//{			
//	struct ComputeTemperatureOp<device::DEVICE_CPU>
//	{
//		void operator()(const rbmd::Id& num_atoms,
//                      const rbmd::Real& mvv2e,
//                      const rbmd::Real* mass,
//                      const rbmd::Real* vx,
//                      const rbmd::Real* vy,
//                      const rbmd::Real* vz,
//                      rbmd::Real& _temp_sum)
//		{
//			for (size_t i = 0; i < num_atoms; i++)
//			{
//				vx += 0.5 * fx / mass_map[i] * dt * fmt2v; 
//				vy += 0.5 * fy / mass_map[i] * dt * fmt2v;
//				vz += 0.5 * fz / mass_map[i] * dt * fmt2v;
//			}
//		}
//	};
//	
//	template struct UpdateTemperatureOp<device::DEVICE_CPU>;
//	template struct UpdateTemperatureOp<device::DEVICE_CPU>;
//}