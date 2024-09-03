#pragma once
#include "common/types.h"
#include "base_reader.h"
#include <memory>

class MmapReader : public BaseReader
{
public:
	MmapReader(const std::string& filePath);
	virtual ~MmapReader();

	int Execute() override;

protected:
	off_t _file_size;
	char* _mapped_memory;
	char* _line_start;
	rbmd::Id _locate;
};