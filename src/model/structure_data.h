#pragma once
#include "structure_header.h"
#include <map>

class StructureData
{
public:
	StructureData()
	{

	}

public:
	StructureHeader _header;
	std::map<rbmd::Id, rbmd::Real> _mass;
	std::map<rbmd::Id, rbmd::Real> _eps;
	std::map<rbmd::Id, rbmd::Real> _sigma;
};