#include "../include/lj_force_field_reader.h"
#include <sstream>
#include <string>
#include "model/md_data.h"
#include "../Utilities/string_util.h"
#include "data_manager.h"
#include <hip/hip_runtime.h>

#define HIP_CHECK(call) { \
    hipError_t err = call; \
    if (err != hipSuccess) { \
        std::cerr << "HIP error: " << hipGetErrorString(err) << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
        exit(err); \
    } \
}
LJForceFieldReader::LJForceFieldReader(MDData& data, const std::string& atom_style, const std::string& force_field, off_t& file_size, char*& mapped_memory, char*& line_start, rbmd::Id& locate) :
	ForceFieldReader(data, atom_style, force_field, file_size, mapped_memory, line_start, locate)
{
}

int LJForceFieldReader::ReadData()
{
	try
	{
		auto& info = _md_data._structure_info_data;
		for (; _locate < _file_size; ++_locate)
		{
			if (_mapped_memory[_locate] == '\n')
			{
				auto line = std::string(_line_start, &_mapped_memory[_locate]);
				std::istringstream iss(line);
				if (rbmd::IsLegalLine(line))
				{
					if (line.find("Masses") != std::string::npos)
					{
						ReadMass(info->_num_atoms_type);
						//std::cout <<"Masses" << std::endl;
					}
					else if (line.find("Pair Coeffs") != std::string::npos)
					{
						ReadPairCoeffs(info->_num_atoms_type);
						//std::cout << "Pair Coeffs" << std::endl;
					}
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

void LJForceFieldReader::AllocateDataSpace()
{
	auto& md_data = DataManager::getInstance().getMDData();
	auto& info = _md_data._structure_info_data;
	auto& data = _md_data._force_field_data;
	HIP_CHECK(hipHostMalloc(&(data->_h_mass), info->_num_atoms_type *sizeof(rbmd::Real)));
	HIP_CHECK(hipHostMalloc(&(data->_h_eps), info->_num_atoms_type * sizeof(rbmd::Real)));
	HIP_CHECK(hipHostMalloc(&(data->_h_sigma), info->_num_atoms_type *sizeof(rbmd::Real)));
}

int LJForceFieldReader::ReadMass(const rbmd::Id& numAtomTypes)
{
	try
	{
		auto force_filed = std::dynamic_pointer_cast<LJForceFieldData>(_md_data._force_field_data);
		auto& mass = force_filed->_h_mass;
		rbmd::Id atom_type;
		rbmd::Real value;

		_line_start = &_mapped_memory[_locate];
		for (auto num = 0; _locate < _file_size && num < numAtomTypes; ++_locate)
		{
			if (_mapped_memory[_locate] == '\n')
			{
				auto line = std::string(_line_start, &_mapped_memory[_locate]);
				std::istringstream iss(line);
				if (rbmd::IsLegalLine(line))
				{
					iss >> atom_type >> value;
					//std::cout << atom_type << " " << value << std::endl;
					mass[atom_type - 1] = value;
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

int LJForceFieldReader::ReadPairCoeffs(const rbmd::Id& numAtomTypes)
{
	try
	{
		auto force_filed = std::dynamic_pointer_cast<LJForceFieldData>(_md_data._force_field_data);
		auto& eps = force_filed->_h_eps;
		auto& sigma = force_filed->_h_sigma;
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
				if (rbmd::IsLegalLine(line))
				{
					iss >> atom_type >> eps_value >> sigma_value;
					//std::cout << atom_type << " " << eps_value << " " << sigma_value << std::endl;
					eps[atom_type - 1] = eps_value;
					sigma[atom_type - 1] = sigma_value;
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
