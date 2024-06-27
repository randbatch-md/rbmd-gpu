#include "md_system.h"
#include "base/device_types.h"
#include "near_force/direct_truncation/direct_truncation_op.h"
#include "base/memory/memory_op.h"

int MDSystem::Evolve()
{
	auto& structure_info = _md_data._structure_info;
	auto& potential_data = _md_data._potential_data;
	auto& structure_data = _md_data._structure_data;
	auto nAtoms = structure_info._num_atoms;
	rbmd::Real* d_dt = nullptr, * d_fmt2v = nullptr, d_mass = nullptr, d_v = nullptr, d_force = nullptr;
	
	//rbmd::Real dt = 0.5, fmt2v = 1.0;
	//rbmd::Real3* force = new rbmd::Real3(nAtoms);
	//op::resize_memory_op(d_dt, 1);
	//op::resize_memory_op(d_fmt2v, 1);
	//op::resize_memory_op(d_mass, potential_data._mass.size());
	//op::resize_memory_op(d_v, nAtoms);
	//op::resize_memory_op(d_force, nAtoms);

	//op::sync_memory_h2d_op<rbmd::Real, device::DEVICE_GPU>(&dt, &d_dt, 1);
	//op::sync_memory_h2d_op<rbmd::Real, device::DEVICE_GPU>(&fmt2v, &d_fmt2v, 1);
	//op::sync_memory_h2d_op<rbmd::Real, device::DEVICE_GPU>(potential_data._mass.data(), &d_mass, 1);
	//op::sync_memory_h2d_op<rbmd::Real3, device::DEVICE_GPU>(structure_data._velocities.data(), &d_v, nAtoms);
	//op::sync_memory_h2d_op<rbmd::Real3, device::DEVICE_GPU>(force, &d_force, nAtoms);

	op::direct_truncation_op<rbmd::Real, device::DEVICE_GPU>()(
		1,
		structure_info._num_atoms,
		d_dt,
		d_fmt2v,
		d_mass,
		d_v,
		d_force);

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
