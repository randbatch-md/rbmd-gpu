#pragma once

#include "mmap_reader.h"
#include "potential_data.h"

class StructureInfo;

class PotentialReader : public MmapReader
{
public:
	PotentialReader(const std::string& filePath, const StructureInfo& info, PotentialData& potential);
	~PotentialReader();

private:
	const StructureInfo& _structure_info;
	PotentialData& _potential_data;
};