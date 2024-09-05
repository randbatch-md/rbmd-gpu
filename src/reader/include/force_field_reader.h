#pragma once

#include "mmap_reader.h"
#include "force_field/force_field_data.h"

class MDData;
class ForceFieldReader : public MmapReader
{
public:
	ForceFieldReader(MDData& data, const std::string& atom_style, const std::string& force_field, off_t& file_size, char*& mapped_memory, char*& line_start, rbmd::Id& locate);
	~ForceFieldReader() = default;

	int Execute() override;
protected:
	virtual int ReadData() = 0;
	virtual void AllocateDataSpace() = 0;

	MDData& _md_data;
};