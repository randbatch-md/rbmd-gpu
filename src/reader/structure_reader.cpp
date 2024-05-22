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
	try
	{
		auto& info = _md_data._structure_info;
		std::vector<int> a;

		for (; _locate < _file_size; ++_locate)
		{
			if (_mapped_memory[_locate] == '\n')
			{
				auto line = std::string(_line_start, &_mapped_memory[_locate]);
				std::istringstream iss(line);
				if (!line.empty() && line != "\n")
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
	}
	catch (const std::exception& e)
	{
		//log
		return -1;
	}

	return 0;
}

int StructureReder::ReadPotential()
{
	try
	{
		auto& info = _md_data._structure_info;
		for (; _locate < _file_size; ++_locate)
		{
			if (_mapped_memory[_locate] == '\n')
			{
				auto line = std::string(_line_start, &_mapped_memory[_locate]);
				std::istringstream iss(line);

				if (!line.empty() && line != "\n")
				{
					if (line.find("Masses") != std::string::npos)
					{
						ReadMass(info._num_atoms_type);
					}
					else if (line.find("Pair Coeffs") != std::string::npos)
					{
						ReadPairCoeffs(info._num_atoms_type);
					}
					else if (line.find("Bond Coeffs") != std::string::npos)
					{
						ReadBondCoeffs(info._num_bound_type);
					}
					else if (line.find("Angle Coeffs") != std::string::npos)
					{
						ReadAngleCoeffs(info._num_angle_type);
					}
					else if (line.find("group") != std::string::npos)
					{
						//reserve
					}
					else if (line.find("Atoms") != std::string::npos)
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

int StructureReder::ReadMass(const rbmd::Id& numAtomTypes)
{
	try
	{
		auto& mass = _md_data._potential_data._mass;
		rbmd::Id atom_type;
		rbmd::Real value;

		_line_start = &_mapped_memory[_locate];
		for (auto num = 0; _locate < _file_size && num < numAtomTypes; ++_locate)
		{
			if (_mapped_memory[_locate] == '\n')
			{
				auto line = std::string(_line_start, &_mapped_memory[_locate]);
				std::istringstream iss(line);
				if (!line.empty() && line != "\n")
				{
					iss >> atom_type >> value;
					mass.insert(std::make_pair(atom_type, value));
					++num;
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

int StructureReder::ReadPairCoeffs(const rbmd::Id& numAtomTypes)
{
	try
	{
		auto& eps = _md_data._potential_data._eps;
		auto& sigma = _md_data._potential_data._sigma;
		rbmd::Id atom_type;
		rbmd::Real eps_value;
		rbmd::Real sigma_value;

		_line_start = &_mapped_memory[_locate];
		for (auto num = 0; _locate < _file_size && num < numAtomTypes; ++_locate)
		{
			if (_mapped_memory[_locate] == '\n')
			{
				auto line = std::string(_line_start, &_mapped_memory[_locate]);
				std::istringstream iss(line);
				if (!line.empty() && line != "\n")
				{
					iss >> atom_type >> eps_value >> sigma_value;
					eps.insert(std::make_pair(atom_type, eps_value));
					sigma.insert(std::make_pair(atom_type, sigma_value));
					++num;
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

int StructureReder::ReadBondCoeffs(const rbmd::Id& numBondTypes)
{
	try
	{

	}
	catch (const std::exception& e)
	{
		//log
		return -1;
	}

	return 0;
}

int StructureReder::ReadAngleCoeffs(const rbmd::Id& numAngleTypes)
{
	try
	{

	}
	catch (const std::exception& e)
	{
		//log
		return -1;
	}

	return 0;
}
