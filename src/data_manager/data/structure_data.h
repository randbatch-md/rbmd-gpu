#pragma once
#include "../model/types.h"
#include "../model/SOA3DArray.h"

struct StructureData
{
	SOA3DArray<rbmd::Real> _h_position;
	std::vector<rbmd::Id> _h_atoms_id;
	std::vector<rbmd::Id> _h_atoms_type;

	SOA3DArray<rbmd::Real> _h_velocity;

	SOA3DArray<rbmd::Real> _h_force;
};
