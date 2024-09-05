#include "../include/cvff_force_field_reader.h"
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
CVFFForceFieldReader::CVFFForceFieldReader(MDData& data, const std::string& atom_style, const std::string& force_field, off_t& file_size, char*& mapped_memory, char*& line_start, rbmd::Id& locate) :
	ForceFieldReader(data, atom_style, force_field, file_size, mapped_memory, line_start, locate)
{
}

int CVFFForceFieldReader::ReadData()
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
					else if (line.find("Bond Coeffs") != std::string::npos)
					{
						ReadBondCoeffs(info->_num_bounds_type);
						//std::cout << "Bond Coeffs" << std::endl;
					}
					else if (line.find("Angle Coeffs") != std::string::npos)
					{
						ReadAngleCoeffs(info->_num_angles_type);
						//std::cout << "Angle Coeffs" << std::endl;
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

void CVFFForceFieldReader::AllocateDataSpace()
{
	auto& md_data = DataManager::getInstance().getMDData();
	auto& info = _md_data._structure_info_data;
	auto& data = _md_data._force_field_data;
	HIP_CHECK(hipHostMalloc(&(data->_h_mass), info->_num_atoms_type * sizeof(rbmd::Real)));
	HIP_CHECK(hipHostMalloc(&(data->_h_eps), info->_num_atoms_type * sizeof(rbmd::Real)));
	HIP_CHECK(hipHostMalloc(&(data->_h_sigma), info->_num_atoms_type * sizeof(rbmd::Real)));
	
	HIP_CHECK(hipHostMalloc(&(data->_h_bond_coeffs_k), info->_num_bounds_type * sizeof(rbmd::Real)));
	HIP_CHECK(hipHostMalloc(&(data->_h_bond_coeffs_equilibrium), info->_num_bounds_type * sizeof(rbmd::Real)));
	
	HIP_CHECK(hipHostMalloc(&(data->_h_angle_coeffs_k), info->_num_angles_type * sizeof(rbmd::Real)));
	HIP_CHECK(hipHostMalloc(&(data->_h_angle_coeffs_equilibrium), info->_num_angles_type * sizeof(rbmd::Real)));
	
	HIP_CHECK(hipHostMalloc(&(data->_h_dihedral_coeffs_k), info->_num_dihedrals_type * sizeof(rbmd::Real)));
	HIP_CHECK(hipHostMalloc(&(data->_h_dihedral_coeffs_sign), info->_num_dihedrals_type * sizeof(rbmd::Id)));
	HIP_CHECK(hipHostMalloc(&(data->_h_dihedral_coeffs_multiplicity), info->_num_dihedrals_type * sizeof(rbmd::Id)));
}

int CVFFForceFieldReader::ReadMass(const rbmd::Id& numAtomTypes)
{
	try
	{
		auto force_filed = std::dynamic_pointer_cast<CVFFForceFieldData>(_md_data._force_field_data);
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

int CVFFForceFieldReader::ReadPairCoeffs(const rbmd::Id& numAtomTypes)
{
	try
	{
		auto force_filed = std::dynamic_pointer_cast<CVFFForceFieldData>(_md_data._force_field_data);
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

int CVFFForceFieldReader::ReadBondCoeffs(const rbmd::Id& numBondTypes)
{
	try
	{
		auto force_filed = std::dynamic_pointer_cast<CVFFForceFieldData>(_md_data._force_field_data);
		auto& bond_coeffs_k = force_filed->_h_bond_coeffs_k;
		auto& bond_coeffs_equilibrium = force_filed ->_h_bond_coeffs_equilibrium;
		rbmd::Id bound_type;
		rbmd::Real bond_coeffs_k_value;
		rbmd::Real equilibrium_value;

		_line_start = &_mapped_memory[_locate];
		for (auto num = 0; _locate < _file_size && num < numBondTypes; ++_locate)
		{
			if (_mapped_memory[_locate] == '\n')
			{
				auto line = std::string(_line_start, &_mapped_memory[_locate]);
				std::istringstream iss(line);
				if (rbmd::IsLegalLine(line))
				{
					iss >> bound_type >> bond_coeffs_k_value >> equilibrium_value;
					//std::cout << bound_type << " " << bond_coeffs_k_value << " " << equilibrium_value << std::endl;
					bond_coeffs_k[bound_type - 1] = bond_coeffs_k_value;
					bond_coeffs_equilibrium[bound_type - 1] = equilibrium_value;
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

int CVFFForceFieldReader::ReadAngleCoeffs(const rbmd::Id& numAngleTypes)
{
	try
	{
		auto force_filed = std::dynamic_pointer_cast<CVFFForceFieldData>(_md_data._force_field_data);
		auto& angle_coeffs_k = force_filed->_h_angle_coeffs_k;
		auto& angle_coeffs_equilibrium = force_filed->_h_angle_coeffs_equilibrium;
		rbmd::Id angle_type;
		rbmd::Real angle_coeffs_k_value;
		rbmd::Real equilibrium_value;

		_line_start = &_mapped_memory[_locate];
		for (auto num = 0; _locate < _file_size && num < numAngleTypes; ++_locate)
		{
			if (_mapped_memory[_locate] == '\n')
			{
				auto line = std::string(_line_start, &_mapped_memory[_locate]);
				std::istringstream iss(line);
				if (rbmd::IsLegalLine(line))
				{
					iss >> angle_type >> angle_coeffs_k_value >> equilibrium_value;
					//std::cout << angle_type << " " << angle_coeffs_k_value << " " << equilibrium_value << std::endl;
					angle_coeffs_k[angle_type - 1] = angle_coeffs_k_value;
					angle_coeffs_equilibrium[angle_type - 1] = equilibrium_value;
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

int CVFFForceFieldReader::ReadDihedralsCoeffs(const rbmd::Id& numAngleTypes)
{
	try
	{
		auto force_filed = std::dynamic_pointer_cast<CVFFForceFieldData>(_md_data._force_field_data);
		auto& dihedral_coeffs_k = force_filed->_h_dihedral_coeffs_k;
		auto& dihedral_coeffs_sign = force_filed->_h_dihedral_coeffs_sign;
		auto& dihedral_coeffs_multiplicity = force_filed->_h_dihedral_coeffs_multiplicity;
		rbmd::Id dihedral_type;
		rbmd::Real dihedral_coeffs_k_value;
		rbmd::Real sign_value;
		rbmd::Real multiplicity_value;

		_line_start = &_mapped_memory[_locate];
		for (auto num = 0; _locate < _file_size && num < numAngleTypes; ++_locate)
		{
			if (_mapped_memory[_locate] == '\n')
			{
				auto line = std::string(_line_start, &_mapped_memory[_locate]);
				std::istringstream iss(line);
				if (rbmd::IsLegalLine(line))
				{
					iss >> dihedral_type >> dihedral_coeffs_k_value >> sign_value >> multiplicity_value;
					//std::cout << angle_type << " " << angle_coeffs_k_value << " " << equilibrium_value << std::endl;
					dihedral_coeffs_k[dihedral_type - 1] = dihedral_coeffs_k_value;
					dihedral_coeffs_sign[dihedral_type - 1] = sign_value;
					dihedral_coeffs_multiplicity[dihedral_type - 1] = multiplicity_value;
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