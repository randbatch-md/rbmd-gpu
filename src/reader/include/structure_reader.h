#pragma once
#include "mmap_reader.h"

struct MDData;
class StructureReder : public MmapReader
{
public:
	StructureReder(const std::string& filePath, MDData& data);
	~StructureReder() = default;

	int Execute() override;

protected:
	virtual int ReadData() = 0;
	virtual void AllocateDataSpace() = 0;

private:
	int ReadHeader();
	int ReadForceField();
	int ReadMass(const rbmd::Id& numAtomTypes);
	int ReadPairCoeffs(const rbmd::Id& numAtomTypes);
	int ReadBondCoeffs(const rbmd::Id& numBondTypes);
	int ReadAngleCoeffs(const rbmd::Id& numAngleTypes);
	//int ReadGroup();

protected:
	MDData& _md_data;
	//std::shared_ptr<ForceFieldReader> _force_field_reader;
};
