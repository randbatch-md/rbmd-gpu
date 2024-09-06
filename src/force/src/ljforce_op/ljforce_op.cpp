//#include "ljforce_op.h"
//#include "device_types.h"
//#include<iostream>
//
//namespace op
//{
//	template <typename FPTYPE>
//	struct ComputeLJforceOp<FPTYPE, device::DEVICE_CPU>
//	{
//		//to do openmp
		//void operator()(Box& box,
		//	            const rbmd::Id& N,
		//	            const rbmd::Id* atoms_type,
		//	            const rbmd::Id* molecular_type,
		//	            const rbmd::Real* sigma,
		//	            const rbmd::Real* eps,
		//	            const rbmd::Id* start_id,
		//                const rbmd::Id* end_id,
		//                const rbmd::Id* id_verletlist,
		//	            const rbmd::Real* px,
		//	            const rbmd::Real* py,
		//	            const rbmd::Real* pz,
		//	            rbmd::Real* force_x,
		//	            rbmd::Real* force_y,
		//	            rbmd::Real* force_z,
		//	            rbmd::Real* evdwl);
//		{
//		}
//	};
//}