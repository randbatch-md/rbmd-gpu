#pragma once

#include <vector>

struct StructureData
{
	std::vector<rbmd::Id> _atom_ids;
	std::vector<rbmd::Id> _atom_types;
	std::vector<rbmd::vec3> _positions;
	std::vector<rbmd::vec3> _velocities;
};