#include "atomic_reader.h"
#include <sstream>
#include <string>

AtomicReader::AtomicReader(const std::string& filePath, MDData& data) :
	StructureReder(filePath, data)
{

}

int AtomicReader::ReadData()
{
	for (; _locate < _file_size; ++_locate)
	{
		if (_mapped_memory[_locate] == '\n')
		{
			auto line = std::string(_line_start, &_mapped_memory[_locate]);
			std::istringstream iss(line);


			_line_start = &_mapped_memory[_locate];
		}
	}

	return 0;
}
