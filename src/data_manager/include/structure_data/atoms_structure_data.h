#pragma once
#include "structure_data.h"

class AtomsStructureData : public StructureData 
{
public:
	bool checkStructure() const override
	{
		return true;
	}
};
