#pragma once
#include "structure_data.h"

class AtomsStructureData : public StructureData 
{
	bool checkStructure() const override
	{
		return true;
	}
};
