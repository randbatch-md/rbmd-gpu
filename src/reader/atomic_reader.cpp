#include "atomic_reader.h"
#include <sstream>
#include <string>
#include "md_data.h"
#include "string_util.h"

AtomicReader::AtomicReader(const std::string& filePath, MDData& data) :
	StructureReder(filePath, data)
{

}

int AtomicReader::ReadData()
{
	try
	{
		auto& num_atoms = _md_data._structure_info._num_atoms;
		//read position
		for (; _locate < _file_size; ++_locate)
		{
			if (_mapped_memory[_locate] == '\n')
			{
				auto line = std::string(_line_start, &_mapped_memory[_locate]);
				std::istringstream iss(line);

				if (line.find("Atoms") != std::string::npos)
				{
					ReadAtoms(num_atoms);
				}
				else if (line.find("Velocities") != std::string::npos)
				{
					ReadVelocity(num_atoms);
				}
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
	data._atom_ids.resize(info._num_atoms);
	data._atom_types.resize(info._num_atoms);
	data._positions.resize(info._num_atoms);
	data._velocities.resize(info._num_atoms);
}

int AtomicReader::ReadAtoms(const rbmd::Id& atoms_num)
{
	try
	{
		auto& ids = _md_data._structure_data._atom_ids;
		auto& types = _md_data._structure_data._atom_types;
		auto& positions = _md_data._structure_data._positions;
		rbmd::Id atom_id;

		_line_start = &_mapped_memory[_locate];
		for (auto num = 0; _locate < _file_size && num < atoms_num; ++_locate)
		{
			if (_mapped_memory[_locate] == '\n')
			{
				auto line = std::string(_line_start, &_mapped_memory[_locate]);
				std::istringstream iss(line);
				if (rbmd::IsLegalLine(line))
				{
					iss >> atom_id;
					auto index = atom_id - 1;
					ids[index] = atom_id;
					iss >> types[index];
					iss >> positions[index].data[0] >> positions[index].data[1] >> positions[index].data[2];
					++num;
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

int AtomicReader::ReadVelocity(const rbmd::Id& atoms_num)
{
	try
	{
		auto& velocities = _md_data._structure_data._velocities;
		rbmd::Id atom_id;

		_line_start = &_mapped_memory[_locate];
		for (auto num = 0; _locate < _file_size && num < atoms_num; ++_locate)
		{
			if (_mapped_memory[_locate] == '\n')
			{
				auto line = std::string(_line_start, &_mapped_memory[_locate]);
				std::istringstream iss(line);
				if (rbmd::IsLegalLine(line))
				{
					iss >> atom_id;
					auto index = atom_id - 1;
					iss >> velocities[index].data[0] >> velocities[index].data[1] >> velocities[index].data[2];
					++num;
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
