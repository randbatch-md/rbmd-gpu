#pragma once
#include "structure_data.h"

struct LJStructureData : StructureData 
{
	SOA3DArrayH<rbmd::Real> _h_velocity;
};
