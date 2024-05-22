#include "atomic_reader.h"
#include <sstream>
#include <string>
#include "md_data.h"

AtomicReader::AtomicReader(const std::string& filePath, MDData& data) :
	StructureReder(filePath, data)
{

}

int AtomicReader::ReadData()
{
	try
	{
		for (; _locate < _file_size; ++_locate)
		{
			if (_mapped_memory[_locate] == '\n')
			{
				auto line = std::string(_line_start, &_mapped_memory[_locate]);
				std::istringstream iss(line);

				if (!line.empty() && line != "\n")
				{
					if (line.find("Atoms") != std::string::npos)
					{

						break;
					}
				}

				_line_start = &_mapped_memory[_locate];
			}
		}

	}
	catch (const std::exception& e)
	{
		//log
		return -1;
	}

	return 0;
}

void AtomicReader::AllocateDataSpace()
{
	auto& info = _md_data._structure_info;
	auto& data = _md_data._structure_data;
	data._position.resize(info._num_atoms);

}
