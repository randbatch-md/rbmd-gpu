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
	int ReadPotential();

protected:
	MDData& _md_data;
};
