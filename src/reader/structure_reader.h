#pragma once
#include "mmap_reader.h"
#include "model/structure_header.h"
#include "model/structure_data.h"

class StructureReder : public MmapReader
{
public:
	StructureReder(const std::string& filePath);
	~StructureReder() = default;

	int Execute() override;

public:
	auto& GetData() { return _Data; }

protected:
	virtual int ReadData() = 0;

private:
	int ReadHeader();

private:
	std::unique_ptr<StructureHeader> _header;
	std::unique_ptr<StructureData> _Data;
};
