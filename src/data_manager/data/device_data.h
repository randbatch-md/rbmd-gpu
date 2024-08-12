#pragma once
#include "../model/soa_3d_array_d.h"
#include "../model/types.h"

struct DeviceData
{
	SOA3DArrayD<rbmd::Real> _h_position;
	thrust::device_vector<rbmd::Id> _h_atoms_id;
	thrust::device_vector<rbmd::Id> _h_atoms_type;

	SOA3DArrayD<rbmd::Real> _h_velocity;

	SOA3DArrayD<rbmd::Real> _h_force;
};
