#pragma once

#include "force_field_data.h"

class CVFFForceFieldData : public ForceFieldData
{
public:
	bool checkForceField() const override
	{
		return true;
	}

	///bond
	rbmd::Real* _h_bond_coeffs_k;
	rbmd::Real* _h_bond_coeffs_equilibrium;

	///angle
	rbmd::Real* _h_angle_coeffs_k;
	rbmd::Real* _h_angle_coeffs_equilibrium;

	///dihedral
	rbmd::Real* _h_dihedral_coeffs_k;
	rbmd::Id* _h_dihedral_coeffs_sign;
	rbmd::Id* _h_dihedral_coeffs_multiplicity;
};
