//#include "../include/charge_reader.h"
//#include <sstream>
//#include <string>
//#include "model/md_data.h"
//#include "../Utilities/string_util.h"
//#include "data_manager.h"
//#include <hip/hip_runtime.h>
//
//#define HIP_CHECK(call) { \
//    hipError_t err = call; \
//    if (err != hipSuccess) { \
//        std::cerr << "HIP error: " << hipGetErrorString(err) << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
//        exit(err); \
//    } \
//}
//Charge_Reader ::Charge_Reader (const std::string& filePath, MDData& data) :
//	StructureReder(filePath, data)
//{
//
//}
//
//int Charge_Reader ::ReadData()
//{
//	try
//	{
//		auto& num_atoms = _md_data._structure_info_data->_num_atoms;
//		read position
//		for (; _locate < _file_size; ++_locate)
//		{
//			if (_mapped_memory[_locate] == '\n')
//			{
//				auto line = std::string(_line_start, &_mapped_memory[_locate]);
//				std::istringstream iss(line);
//
//				if (line.find("Atoms") != std::string::npos)
//				{
//					ReadAtoms(num_atoms);
//					std::cout << "Atoms" << std::endl;
//				}
//				else if (line.find("Velocities") != std::string::npos)
//				{
//					ReadVelocity(num_atoms);
//				}
//			}
//		}
//
//	}
//	catch (const std::exception& e)
//	{
//		log
//		return -1;
//	}
//
//	return 0;
//}
//
//void Charge_Reader ::AllocateDataSpace()
//{
//	auto& md_data = DataManager::getInstance().getMDData();
//	auto& info = _md_data._structure_info_data;
//	auto& data = _md_data._structure_data;
//	HIP_CHECK(hipHostMalloc(&(data->_h_atoms_id), info->_num_atoms*sizeof(rbmd::Id)));
//	HIP_CHECK(hipHostMalloc(&(data->_h_atoms_type), info->_num_atoms * sizeof(rbmd::Id)));
//	HIP_CHECK(hipHostMalloc(&(data->_h_px), info->_num_atoms*sizeof(rbmd::Id)));
//	HIP_CHECK(hipHostMalloc(&(data->_h_py), info->_num_atoms*sizeof(rbmd::Id)));
//	HIP_CHECK(hipHostMalloc(&(data->_h_pz), info->_num_atoms*sizeof(rbmd::Id)));
//	HIP_CHECK(hipHostMalloc(&(data->_h_vx), info->_num_atoms*sizeof(rbmd::Id)));
//	HIP_CHECK(hipHostMalloc(&(data->_h_vy), info->_num_atoms*sizeof(rbmd::Id)));
//	HIP_CHECK(hipHostMalloc(&(data->_h_vz), info->_num_atoms*sizeof(rbmd::Id)));
//}
//
//int Charge_Reader ::ReadAtoms(const rbmd::Id& atoms_num)
//{
//	try
//	{
//		auto& ids = _md_data._structure_data->_h_atoms_id;
//		auto& types = _md_data._structure_data->_h_atoms_type;
//		rbmd::Id atom_id;
//
//		_line_start = &_mapped_memory[_locate];
//		for (auto num = 0; _locate < _file_size && num < atoms_num; ++_locate)
//		{
//			if (_mapped_memory[_locate] == '\n')
//			{
//				auto line = std::string(_line_start, &_mapped_memory[_locate]);
//				std::istringstream iss(line);
//				if (rbmd::IsLegalLine(line))
//				{
//					iss >> atom_id;
//					auto index = atom_id - 1;
//					ids[index] = atom_id;
//					iss >> types[index];
//					iss >> _md_data._structure_data->_h_px[index] >> _md_data._structure_data->_h_py[index] >>_md_data._structure_data->_h_pz[index];
//					++num;
//					std::cout << atom_id << " " << types[index] << " " << _md_data._structure_data->_h_px[index] << " " << _md_data._structure_data->_h_py[index] << " " << _md_data._structure_data->_h_pz[index] << std::endl;
//				}
//				_line_start = &_mapped_memory[_locate];
//			}
//		}
//	}
//	catch (const std::exception& e)
//	{
//		log
//		return -1;
//	}
//
//	return 0;
//}
//
//int Charge_Reader ::ReadVelocity(const rbmd::Id& atoms_num)
//{
//	try
//	{
//		rbmd::Id atom_id;
//
//		_line_start = &_mapped_memory[_locate];
//		for (auto num = 0; _locate < _file_size && num < atoms_num; ++_locate)
//		{
//			if (_mapped_memory[_locate] == '\n')
//			{
//				auto line = std::string(_line_start, &_mapped_memory[_locate]);
//				std::istringstream iss(line);
//				if (rbmd::IsLegalLine(line))
//				{
//					iss >> atom_id;
//					auto index = atom_id - 1;
//					iss >> _md_data._structure_data->_h_vx[index] >> _md_data._structure_data->_h_vy[index] >> _md_data._structure_data->_h_vz[index];
//					++num;
//					std::cout << atom_id << " " << _md_data._structure_data->_h_vx[index] << " " << _md_data._structure_data->_h_vy[index] << " " << _md_data._structure_data->_h_vz[index] << std::endl;
//				}
//				_line_start = &_mapped_memory[_locate];
//			}
//		}
//	}
//	catch (const std::exception& e)
//	{
//		log
//		return -1;
//	}
//
//	return 0;
//}
