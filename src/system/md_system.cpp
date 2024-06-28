#include "md_system.h"
#include "base/device_types.h"
#include "near_force/direct_truncation/direct_truncation_op.h"
#include "base/memory/memory_op.h"

int MDSystem::Evolve()
{
	auto& structure_info = _md_data._structure_info;
	auto& potential_data = _md_data._potential_data;
	auto& structure_data = _md_data._structure_data;
	auto& nAtoms = structure_info._num_atoms;
	rbmd::Real* d_dt = nullptr, * d_fmt2v = nullptr, * d_mass = nullptr;
	rbmd::Real3* d_v = nullptr, * d_force = nullptr;
	
	rbmd::Real dt = 0.5, fmt2v = 1.0;
	rbmd::Real3* force = new rbmd::Real3[nAtoms];
	op::resize_memory_op<rbmd::Real,device::DEVICE_GPU>()(d_dt, 1);
	op::resize_memory_op<rbmd::Real,device::DEVICE_GPU>()(d_fmt2v, 1);
	op::resize_memory_op<rbmd::Real,device::DEVICE_GPU>()(d_mass, potential_data._mass.size());
	op::resize_memory_op<rbmd::Real3,device::DEVICE_GPU>()(d_v, nAtoms);
	op::resize_memory_op<rbmd::Real3,device::DEVICE_GPU>()(d_force, nAtoms);

	op::sync_memory_h2d_op<rbmd::Real, device::DEVICE_GPU>()(d_dt, &dt, 1);
	op::sync_memory_h2d_op<rbmd::Real, device::DEVICE_GPU>()(d_fmt2v, &fmt2v, 1);
	op::sync_memory_h2d_op<rbmd::Real, device::DEVICE_GPU>()(d_mass, potential_data._mass.data(), 1);
	op::sync_memory_h2d_op<rbmd::Real3, device::DEVICE_GPU>()(d_v, structure_data._velocities.data(), nAtoms);
	op::sync_memory_h2d_op<rbmd::Real3, device::DEVICE_GPU>()(d_force, force, nAtoms);

	op::direct_truncation_op<rbmd::Real, device::DEVICE_GPU>()(
		1,
		structure_info._num_atoms,
		d_dt,
		d_fmt2v,
		d_mass,
		d_v,
		d_force);

	//op::delete_memory_op<rbmd::Real, device::DEVICE_GPU>(d_dt);
	//op::delete_memory_op<rbmd::Real, device::DEVICE_GPU>(d_fmt2v);
	//op::delete_memory_op<rbmd::Real, device::DEVICE_GPU>(d_mass);
	//op::delete_memory_op<rbmd::Real3, device::DEVICE_GPU>(d_v);
	//op::delete_memory_op<rbmd::Real3, device::DEVICE_GPU>(d_force);

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
