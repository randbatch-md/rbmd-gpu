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

	rbmd::Id* _h_atoms_id;
	rbmd::Id* _h_atoms_type;
	rbmd::Id* _h_molecular_id;

	///velocity on host
	rbmd::Real* _h_vx;
	rbmd::Real* _h_vy;
	rbmd::Real* _h_vz;

	///force on host
	rbmd::Real* _h_fx;
	rbmd::Real* _h_fy;
	rbmd::Real* _h_fz;
};
