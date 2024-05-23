#pragma once

#include <vector>

struct Atoms
{
	std::vector<rbmd::Id> _ids;
	std::vector<rbmd::Id> _types;
	std::vector<rbmd::vec3> _positions;
};

struct StructureData
{
	Atoms _atoms;
	std::vector<rbmd::vec3> _velocities;
};