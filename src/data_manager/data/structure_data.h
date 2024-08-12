#pragma once
#include "../model/types.h"
#include "../model/soa_3d_array_h.h"
#include "object.h"

class StructureData : public Object
{
public:
	virtual bool checkData()const = 0;

	SOA3DArrayH<rbmd::Real> _h_position;
	std::vector<rbmd::Id> _h_atoms_id;
	std::vector<rbmd::Id> _h_atoms_type;
	std::vector<rbmd::Id> _h_molecular_id;

	SOA3DArrayH<rbmd::Real> _h_velocity;

	SOA3DArrayH<rbmd::Real> _h_force;
};
