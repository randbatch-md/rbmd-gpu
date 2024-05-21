#include "potential_reader.h"

PotentialReader::PotentialReader(const std::string& filePath, const StructureInfo& info, PotentialData& potential) :
	MmapReader(filePath),
	_structure_info(info),
	_potential_data(potential)
{

}

PotentialReader::~PotentialReader()
{

}
