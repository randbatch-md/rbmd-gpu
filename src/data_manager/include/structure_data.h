#pragma once
#include "types.h"
#include "object.h"

class StructureData : public Object
{
public:
	virtual bool checkStructure()const = 0;

	std::vector<rbmd::Real> _h_px;
	std::vector<rbmd::Real> _h_py;
	std::vector<rbmd::Real> _h_pz;

	std::vector<rbmd::Id> _h_atoms_id;
	std::vector<rbmd::Id> _h_atoms_type;
	std::vector<rbmd::Id> _h_molecular_id;

	std::vector<rbmd::Real> _h_vx;
	std::vector<rbmd::Real> _h_vy;
	std::vector<rbmd::Real> _h_vz;

	std::vector<rbmd::Real> _h_fx;
	std::vector<rbmd::Real> _h_fy;
	std::vector<rbmd::Real> _h_fz;
};
