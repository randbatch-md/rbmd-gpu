#pragma once
#include "basic_structure_data.h"

class FullStructureData : public BasicStructureData
{
public:
	bool checkStructure() const override
	{
		if (false == BasicStructureData::checkStructure())
		{
			return false;
		}

		//
		return true;
	}

	///bond
	rbmd::Id* _h_bond_type;
	rbmd::Id* _h_bond_id0;
	rbmd::Id* _h_bond_id1;

	///angle
	rbmd::Id* _h_angle_type;
	rbmd::Id* _h_angle_id0;
	rbmd::Id* _h_angle_id1;
	rbmd::Id* _h_angle_id2;

	///dihedral
	rbmd::Id* _h_dihedral_type;
	rbmd::Id* _h_dihedral_id0;
	rbmd::Id* _h_dihedral_id1;
	rbmd::Id* _h_dihedral_id2;
	rbmd::Id* _h_dihedral_id3;
};