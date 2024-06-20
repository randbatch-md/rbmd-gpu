#include "md_system.h"
#include "base/device_types.h"
#include "near_force/direct_truncation/direct_truncation_op.h"

int MDSystem::Evolve()
{
	op::direct_truncation_op<rbmd::Real, device::DEVICE_GPU>()();
	//op::LJ();
	return 0;
}

int MDSystem::PreSolve()
{
	return 0;
}

int MDSystem::Solve()
{
	return 0;
}

int MDSystem::PostSolve()
{
	return 0;
}
