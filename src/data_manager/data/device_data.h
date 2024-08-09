#pragma once
#include "../model/SOA3DArrayD.h"
#include "../model/types.h"

struct DeviceData
{
	SOA3DArrayD<rbmd::Real> _h_position;
	thrust::device_vector<rbmd::Id> _h_atoms_id;
	thrust::device_vector<rbmd::Id> _h_atoms_type;

	SOA3DArrayD<rbmd::Real> _h_velocity;

	SOA3DArrayD<rbmd::Real> _h_force;
};
