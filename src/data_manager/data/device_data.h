#pragma once
#include "../model/soa_3d_array_d.h"
#include "../model/types.h"

class DeviceData
{
	SOA3DArrayD<rbmd::Real> _d_position;
	thrust::device_vector<rbmd::Id> _d_atoms_id;
	thrust::device_vector<rbmd::Id> _d_atoms_type;
	thrust::device_vector<rbmd::Id> _d_molecular_id;

	SOA3DArrayD<rbmd::Real> _d_velocity;

	SOA3DArrayD<rbmd::Real> _d_force;
};
