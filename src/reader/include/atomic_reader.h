#pragma once
#include "structure_reader.h"
#include "../common/types.h"

class AtomicReader : public StructureReder
{
public:
	AtomicReader(const std::string& filePath, MDData& data);
	~AtomicReader() = default;

protected:
	int ReadData() override;
	void AllocateDataSpace() override;

private:
	int ReadAtoms(const rbmd::Id& atoms_num);
	int ReadVelocity(const rbmd::Id& atoms_num);
};