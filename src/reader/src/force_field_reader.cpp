#include "../include/force_field_reader.h"

ForceFieldReader::ForceFieldReader(MDData& data, const std::string& atom_style, const std::string& force_field, off_t& file_size, char*& mapped_memory, char*& line_start, rbmd::Id& locate) :
	MmapReader(" ", atom_style, force_field),
	_md_data(data)
{
	_file_size = file_size;
	_mapped_memory = mapped_memory;
	_line_start = line_start;
	_locate = locate;
}

int ForceFieldReader::Execute()
{
	AllocateDataSpace();
	if (-1 == ReadData())
	{
		//log
		return -1;
	}
	return 0;
}
