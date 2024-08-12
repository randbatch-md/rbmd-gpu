#pragma once
#include "../model/types.h"
#include "../model/soa_3d_array_h.h"

struct StructureData
{
	SOA3DArrayH<rbmd::Real> _h_position;
	std::vector<rbmd::Id> _h_atoms_id;
	std::vector<rbmd::Id> _h_atoms_type;

	SOA3DArrayH<rbmd::Real> _h_force;
};
