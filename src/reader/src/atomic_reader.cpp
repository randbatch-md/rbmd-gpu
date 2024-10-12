#include "../include/atomic_reader.h"

#include <hip/hip_runtime.h>

#include <sstream>
#include <string>

#include "../Utilities/string_util.h"
#include "data_manager.h"
#include "model/md_data.h"

#define HIP_CHECK(call)                                              \
  {                                                                  \
    hipError_t err = call;                                           \
    if (err != hipSuccess) {                                         \
      std::cerr << "HIP error: " << hipGetErrorString(err) << " at " \
                << __FILE__ << ":" << __LINE__ << std::endl;         \
      exit(err);                                                     \
    }                                                                \
  }
AtomicReader::AtomicReader(const std::string& filePath, MDData& data)
    : StructureReder(filePath, data) {}

int AtomicReader::ReadData() {
  try {
    auto& num_atoms = _md_data._structure_info_data->_num_atoms;
    // read position
    for (; _locate < _file_size; ++_locate) {
      if (_mapped_memory[_locate] == '\n') {
        auto line = std::string(_line_start, &_mapped_memory[_locate]);
        std::istringstream iss(line);

        if (line.find("Atoms") != std::string::npos) {
          // std::cout << "Atoms" << std::endl;
          ReadAtoms(*num_atoms);
        } else if (line.find("Velocities") != std::string::npos) {
          // std::cout << "Velocities" << std::endl;
          ReadVelocity(*num_atoms);
        }
      }
    }

  } catch (const std::exception& e) {
    // log
    return -1;
  }

  return 0;
}

void AtomicReader::AllocateDataSpace() {
  auto atom_style =
      DataManager::getInstance().getConfigData()->Get<std::string>(
          "atom_style", "init_configuration", "read_data");
  if ("atomic" == atom_style) {
    auto& md_data = DataManager::getInstance().getMDData();
    auto& info = _md_data._structure_info_data;
    auto& data = _md_data._structure_data;

    HIP_CHECK(MALLOCHOST(&(data->_h_atoms_id),
                         *(info->_num_atoms) * sizeof(rbmd::Id)));
    HIP_CHECK(MALLOCHOST(&(data->_h_atoms_type),
                         *(info->_num_atoms) * sizeof(rbmd::Id)));
    HIP_CHECK(
        MALLOCHOST(&(data->_h_px), *(info->_num_atoms) * sizeof(rbmd::Id)));
    HIP_CHECK(
        MALLOCHOST(&(data->_h_py), *(info->_num_atoms) * sizeof(rbmd::Id)));
    HIP_CHECK(
        MALLOCHOST(&(data->_h_pz), *(info->_num_atoms) * sizeof(rbmd::Id)));
    HIP_CHECK(
        MALLOCHOST(&(data->_h_vx), *(info->_num_atoms) * sizeof(rbmd::Id)));
    HIP_CHECK(
        MALLOCHOST(&(data->_h_vy), *(info->_num_atoms) * sizeof(rbmd::Id)));
    HIP_CHECK(
        MALLOCHOST(&(data->_h_vz), *(info->_num_atoms) * sizeof(rbmd::Id)));

  } else if ("charge" == atom_style) {
    auto& md_data = DataManager::getInstance().getMDData();
    auto& info = _md_data._structure_info_data;
    auto& charge_structure_data = _md_data._structure_data;
    ChargeStructureData* data =
        dynamic_cast<ChargeStructureData*>(charge_structure_data.get());
    HIP_CHECK(MALLOCHOST(&(data->_h_atoms_id),
                         *(info->_num_atoms) * sizeof(rbmd::Id)));
    HIP_CHECK(MALLOCHOST(&(data->_h_atoms_type),
                         *(info->_num_atoms) * sizeof(rbmd::Id)));
    HIP_CHECK(
        MALLOCHOST(&(data->_h_px), *(info->_num_atoms) * sizeof(rbmd::Id)));
    HIP_CHECK(
        MALLOCHOST(&(data->_h_py), *(info->_num_atoms) * sizeof(rbmd::Id)));
    HIP_CHECK(
        MALLOCHOST(&(data->_h_pz), *(info->_num_atoms) * sizeof(rbmd::Id)));
    HIP_CHECK(
        MALLOCHOST(&(data->_h_vx), *(info->_num_atoms) * sizeof(rbmd::Id)));
    HIP_CHECK(
        MALLOCHOST(&(data->_h_vy), *(info->_num_atoms) * sizeof(rbmd::Id)));
    HIP_CHECK(
        MALLOCHOST(&(data->_h_vz), *(info->_num_atoms) * sizeof(rbmd::Id)));
    HIP_CHECK(
        MALLOCHOST(&(data->_h_charge), *(info->_num_atoms) * sizeof(rbmd::Id)));
  }
}

int AtomicReader::ReadAtoms(const rbmd::Id& atoms_num) {
  try {
    auto& ids = _md_data._structure_data->_h_atoms_id;
    auto& types = _md_data._structure_data->_h_atoms_type;
    rbmd::Id atom_id;
    rbmd::Id atom_type;

    _line_start = &_mapped_memory[_locate];
    auto atom_style =
        DataManager::getInstance().getConfigData()->Get<std::string>(
            "atom_style", "init_configuration", "read_data");
    if ("atomic" == atom_style) {
      for (auto num = 0; _locate < _file_size && num < atoms_num; ++_locate) {
        if (_mapped_memory[_locate] == '\n') {
          auto line = std::string(_line_start, &_mapped_memory[_locate]);
          std::istringstream iss(line);
          if (rbmd::IsLegalLine(line)) {
            iss >> atom_id;
            auto index = atom_id - 1;
            ids[index] = atom_id - 1;
            iss >> atom_type;
            types[index] = atom_type - 1;
            iss >> _md_data._structure_data->_h_px[index] >>
                _md_data._structure_data->_h_py[index] >>
                _md_data._structure_data->_h_pz[index];
            ++num;
            // std::cout << atom_id << " " << types[index] << " " <<
            // _md_data._structure_data->_h_px[index] << " " <<
            // _md_data._structure_data->_h_py[index] << " " <<
            // _md_data._structure_data->_h_pz[index] << std::endl;
          }
          _line_start = &_mapped_memory[_locate];
        }
      }
    } else if ("charge" == atom_style) {
      auto& charge_structure_data = _md_data._structure_data;
      ChargeStructureData* data =
          dynamic_cast<ChargeStructureData*>(charge_structure_data.get());
      for (auto num = 0; _locate < _file_size && num < atoms_num; ++_locate) {
        if (_mapped_memory[_locate] == '\n') {
          auto line = std::string(_line_start, &_mapped_memory[_locate]);
          std::istringstream iss(line);
          if (rbmd::IsLegalLine(line)) {
            iss >> atom_id;
            auto index = atom_id - 1;
            ids[index] = atom_id;
            iss >> types[index] >> data->_h_charge[index];
            iss >> data->_h_px[index] >> data->_h_py[index] >>
                data->_h_pz[index];
            ++num;
            // std::cout << atom_id << " " << types[index] << " " <<
            // data->_h_charge[index]  << " " << data->_h_px[index] << " " <<
            // data->_h_py[index] << " " << data->_h_pz[index] << std::endl;
          }
          _line_start = &_mapped_memory[_locate];
        }
      }
    }
  } catch (const std::exception& e) {
    // log
    return -1;
  }

  return 0;
}

int AtomicReader::ReadVelocity(const rbmd::Id& atoms_num) {
  try {
    rbmd::Id atom_id;

    _line_start = &_mapped_memory[_locate];
    for (auto num = 0; _locate < _file_size && num < atoms_num; ++_locate) {
      if (_mapped_memory[_locate] == '\n') {
        auto line = std::string(_line_start, &_mapped_memory[_locate]);
        std::istringstream iss(line);
        if (rbmd::IsLegalLine(line)) {
          iss >> atom_id;
          auto index = atom_id - 1;
          iss >> _md_data._structure_data->_h_vx[index] >>
              _md_data._structure_data->_h_vy[index] >>
              _md_data._structure_data->_h_vz[index];
          ++num;
          // std::cout << atom_id << " " <<
          // _md_data._structure_data->_h_vx[index] << " " <<
          // _md_data._structure_data->_h_vy[index] << " " <<
          // _md_data._structure_data->_h_vz[index] << std::endl;
        }
        _line_start = &_mapped_memory[_locate];
      }
    }
  } catch (const std::exception& e) {
    // log
    return -1;
  }

  return 0;
}
