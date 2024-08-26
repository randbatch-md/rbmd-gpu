#pragma once

#include "types.h"

class StructureInfoData
{
public:
	StructureInfoData() :
		_num_atoms(0),
		_num_bonds(0),
		_num_angles(0),
		_num_dihedrals(0),
		_num_impropers(0),
		_num_atoms_type(0),
		_num_bounds_type(0),
		_num_angles_type(0)
	{

	}

public:
	rbmd::Id _num_atoms;
	rbmd::Id _num_bonds;
	rbmd::Id _num_angles;
	rbmd::Id _num_dihedrals;
	rbmd::Id _num_impropers;
	rbmd::Id _num_atoms_type;
	rbmd::Id _num_bounds_type;
	rbmd::Id _num_angles_type;
	rbmd::Id _num_dihedrals_type;

	rbmd::Range _range;
};
