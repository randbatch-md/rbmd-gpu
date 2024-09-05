#pragma once

#include "../include/force_field_reader.h"

class LJForceFieldReader : public ForceFieldReader
{
public:
	LJForceFieldReader(MDData& data, const std::string& atom_style, const std::string& force_field, off_t& file_size, char*& mapped_memory, char*& line_start, rbmd::Id& locate);
	~LJForceFieldReader() = default;

protected:
	int ReadData() override;
	void AllocateDataSpace() override;

private:
	int ReadMass(const rbmd::Id& numAtomTypes);
	int ReadPairCoeffs(const rbmd::Id& numAtomTypes);
};