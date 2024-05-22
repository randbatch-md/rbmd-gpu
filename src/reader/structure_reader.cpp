#include "structure_reader.h"
#include <sstream>
#include "md_data.h"

StructureReder::StructureReder(const std::string& filePath, MDData& data) :
	MmapReader(filePath),
	_md_data(data)
{

}

int StructureReder::Execute()
{
	if (-1 == MmapReader::Execute())
	{
		//log
		return -1;
	}

	if (-1 == ReadHeader())
	{
		//log
		return -1;
	}

	if (-1 == ReadPotential())
	{
		//log
		return -1;
	}

	AllocateDataSpace();
	if (-1 == ReadData())
	{
		//log
		return -1;
	}


	return 0;
}

int StructureReder::ReadHeader()
{
	auto& info = _md_data._structure_info;
	std::vector<int> a;

	for (; _locate < _file_size; ++_locate)
	{
		if (_mapped_memory[_locate] == '\n')
		{
			auto line = std::string(_line_start, &_mapped_memory[_locate]);
			std::istringstream iss(line);
			if (!line.empty() && line != "\r")
			{
				if (line.find("atoms") != std::string::npos)
				{
					iss >> info._num_atoms;
				}
				else if (line.find("bonds") != std::string::npos)
				{
					iss >> info._num_bonds;
				}
				else if (line.find("angles") != std::string::npos)
				{
					iss >> info._num_angles;
				}
				else if (line.find("dihedrals") != std::string::npos)
				{
					iss >> info._num_dihedrals;
				}
				else if (line.find("impropers") != std::string::npos)
				{
					iss >> info._num_impropers;
				}
				else if (line.find("atom types") != std::string::npos)
				{
					iss >> info._num_atoms_type;
				}
				else if (line.find("bond types") != std::string::npos)
				{
					iss >> info._num_bound_type;

				}
				else if (line.find("angle types") != std::string::npos)
				{
					iss >> info._num_angle_type;

				}
				else if (line.find("xlo xhi") != std::string::npos)
				{
					iss >> info._range[0][0] >> info._range[0][1];

				}
				else if (line.find("ylo yhi") != std::string::npos)
				{
					iss >> info._range[1][0] >> info._range[1][1];

				}
				else if (line.find("zlo zhi") != std::string::npos)
				{
					iss >> info._range[2][0] >> info._range[2][1];
					_line_start = &_mapped_memory[_locate];
					break;
				}
			}

			_line_start = &_mapped_memory[_locate];
		}
	}

	return 0;
}

int StructureReder::ReadPotential()
{
	for (; _locate < _file_size; ++_locate)
	{
		if (_mapped_memory[_locate] == '\n')
		{
			auto line = std::string(_line_start, &_mapped_memory[_locate]);
			std::istringstream iss(line);

			if (!line.empty() && line != "\r")
			{
				if (line.find("Masses") != std::string::npos)
				{

				}
				else if (line.find("Pair Coeffs") != std::string::npos)
				{

				}
				else if (line.find("Bond Coeffs") != std::string::npos)
				{

				}
				else if (line.find("Angle Coeffs") != std::string::npos)
				{

				}
				else if (line.find("group") != std::string::npos)
				{

				}
				else if (line.find("Atoms") != std::string::npos)
				{
					break;
				}
			}

			_line_start = &_mapped_memory[_locate];
		}
	}

	return 0;
}
