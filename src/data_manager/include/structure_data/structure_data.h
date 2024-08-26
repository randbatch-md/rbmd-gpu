#pragma once
#include "types.h"
#include "object.h"

class StructureData : public Object
{
public:
	virtual bool checkStructure()const = 0;

	///position on host
	rbmd::Real* _h_px;
	rbmd::Real* _h_py;
	rbmd::Real* _h_pz;

	///atoms is on host
	rbmd::Id* _h_atoms_id;

	///atoms type on host
	rbmd::Id* _h_atoms_type;

	///molecular id on host
	rbmd::Id* _h_molecular_id;

	//atoms flag on host
	rbmd::Id* _h_flagX;
	rbmd::Id* _h_flagY;
	rbmd::Id* _h_flagZ;

	///velocity on host
	rbmd::Real* _h_vx;
	rbmd::Real* _h_vy;
	rbmd::Real* _h_vz;

	///force on host
	rbmd::Real* _h_fx;
	rbmd::Real* _h_fy;
	rbmd::Real* _h_fz;
};
