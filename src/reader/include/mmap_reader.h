#pragma once
#include "types.h"
#include "base_reader.h"
#include <memory>

class MmapReader : public BaseReader
{
public:
	MmapReader(const std::string& filePath, const std::string& atom_style, const std::string& force_field);
	virtual ~MmapReader();

	int Execute() override;

protected:
	off_t _file_size;
	char* _mapped_memory;
	char* _line_start;
	rbmd::Id _locate;

	std::string _atom_style;
	std::string _force_field;
};