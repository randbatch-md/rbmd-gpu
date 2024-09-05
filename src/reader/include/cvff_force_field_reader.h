#pragma once

#include "../include/force_field_reader.h"

class CVFFForceFieldReader : public ForceFieldReader
{
public:
	CVFFForceFieldReader(MDData& data, const std::string& atom_style, const std::string& force_field, off_t& file_size, char*& mapped_memory, char*& line_start, rbmd::Id& locate);
	~CVFFForceFieldReader() = default;

protected:
	int ReadData() override;
	void AllocateDataSpace() override;

private:
	int ReadMass(const rbmd::Id& numAtomTypes);
	int ReadPairCoeffs(const rbmd::Id& numAtomTypes);
	int ReadBondCoeffs(const rbmd::Id& numBondTypes);
	int ReadAngleCoeffs(const rbmd::Id& numAngleTypes);
	int ReadDihedralsCoeffs(const rbmd::Id& numAngleTypes);
};