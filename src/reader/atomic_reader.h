#pragma once
#include "structure_reader.h"

class AtomicReader : public StructureReder
{
public:
	AtomicReader(const std::string& filePath, MDData& data);
	~AtomicReader() = default;

protected:
	int ReadData() override;
};