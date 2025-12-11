#include "../include/atomic_reader.h"

#include <sstream>
#include <string>
#include <unordered_set>

#include "../Utilities/string_util.h"
#include "data_manager.h"
#include "model/md_data.h"
#include "output/include/Logger.hpp"
#include "rbmd_define.h"

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
        } else if (line.find("Bonds") != std::string::npos) {
          // std::cout << "Bonds" << std::endl;
          ReadBond(*(_md_data._structure_info_data->_num_bonds));
        } else if (line.find("Angles") != std::string::npos) {
          // std::cout << "Angles" << std::endl;
          ReadAngle(*(_md_data._structure_info_data->_num_angles));
        } else if (line.find("Dihedrals") != std::string::npos) {
          // std::cout << "Dihedrals" << std::endl;
          ReadDihedrals(*(_md_data._structure_info_data->_num_dihedrals));
        } else if (line.find("Impropers") != std::string::npos) {
          //std::cout << "Dihedrals" << std::endl;
          ReadImpropers(*(_md_data._structure_info_data->_num_impropers));
          //std::cout << "Impropers" << std::endl;
        } else if (line.find("Velocities") != std::string::npos) {
          // std::cout << "Velocities" << std::endl;
          ReadVelocity(*num_atoms);
        }
      }
    }
    //SetSpecialBonds_fix();

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
  auto& md_data = DataManager::getInstance().getMDData();
  auto& info = _md_data._structure_info_data;
  if ("atomic" == atom_style) {
    auto& data = _md_data._structure_data;

    CHECK_RUNTIME(MALLOCHOST(&(data->_h_atoms_id),*(info->_num_atoms) * sizeof(rbmd::Id)));
    CHECK_RUNTIME(MALLOCHOST(&(data->_h_atoms_type),
                         *(info->_num_atoms) * sizeof(rbmd::Id)));
    CHECK_RUNTIME(
        MALLOCHOST(&(data->_h_px), *(info->_num_atoms) * sizeof(rbmd::Real)));
    CHECK_RUNTIME(
        MALLOCHOST(&(data->_h_py), *(info->_num_atoms) * sizeof(rbmd::Real)));
    CHECK_RUNTIME(
        MALLOCHOST(&(data->_h_pz), *(info->_num_atoms) * sizeof(rbmd::Real)));
    CHECK_RUNTIME(
        MALLOCHOST(&(data->_h_vx), *(info->_num_atoms) * sizeof(rbmd::Real)));
    CHECK_RUNTIME(
        MALLOCHOST(&(data->_h_vy), *(info->_num_atoms) * sizeof(rbmd::Real)));
    CHECK_RUNTIME(
        MALLOCHOST(&(data->_h_vz), *(info->_num_atoms) * sizeof(rbmd::Real)));

    CHECK_RUNTIME(
    MALLOCHOST(&(data->_h_flagX), *(info->_num_atoms) * sizeof(rbmd::Id)));
    CHECK_RUNTIME(
        MALLOCHOST(&(data->_h_flagY), *(info->_num_atoms) * sizeof(rbmd::Id)));
    CHECK_RUNTIME(
      MALLOCHOST(&(data->_h_flagZ), *(info->_num_atoms) * sizeof(rbmd::Id)));

  } else if ("charge" == atom_style) {
    auto& charge_structure_data = _md_data._structure_data;
    ChargeStructureData* data =
        dynamic_cast<ChargeStructureData*>(charge_structure_data.get());
    CHECK_RUNTIME(MALLOCHOST(&(data->_h_atoms_id),
                         *(info->_num_atoms) * sizeof(rbmd::Id)));
    CHECK_RUNTIME(MALLOCHOST(&(data->_h_atoms_type),
                         *(info->_num_atoms) * sizeof(rbmd::Id)));
    CHECK_RUNTIME(
        MALLOCHOST(&(data->_h_px), *(info->_num_atoms) * sizeof(rbmd::Real)));
    CHECK_RUNTIME(
        MALLOCHOST(&(data->_h_py), *(info->_num_atoms) * sizeof(rbmd::Real)));
    CHECK_RUNTIME(
        MALLOCHOST(&(data->_h_pz), *(info->_num_atoms) * sizeof(rbmd::Real)));
    CHECK_RUNTIME(
        MALLOCHOST(&(data->_h_vx), *(info->_num_atoms) * sizeof(rbmd::Real)));
    CHECK_RUNTIME(
        MALLOCHOST(&(data->_h_vy), *(info->_num_atoms) * sizeof(rbmd::Real)));
    CHECK_RUNTIME(
        MALLOCHOST(&(data->_h_vz), *(info->_num_atoms) * sizeof(rbmd::Real)));
    CHECK_RUNTIME(
        MALLOCHOST(&(data->_h_charge), *(info->_num_atoms) * sizeof(rbmd::Real)));

    CHECK_RUNTIME(
    MALLOCHOST(&(data->_h_flagX), *(info->_num_atoms) * sizeof(rbmd::Id)));
    CHECK_RUNTIME(
        MALLOCHOST(&(data->_h_flagY), *(info->_num_atoms) * sizeof(rbmd::Id)));
    CHECK_RUNTIME(
      MALLOCHOST(&(data->_h_flagZ), *(info->_num_atoms) * sizeof(rbmd::Id)));

  } else if ("full" == atom_style) {
      auto& full_structure_data = _md_data._structure_data;
      FullStructureData* data =
          dynamic_cast<FullStructureData*>(full_structure_data.get());
      CHECK_RUNTIME(MALLOCHOST(&(data->_h_atoms_id),
          *(info->_num_atoms) * sizeof(rbmd::Id)));
      CHECK_RUNTIME(MALLOCHOST(&(data->_h_atoms_type),
          *(info->_num_atoms) * sizeof(rbmd::Id)));
      CHECK_RUNTIME(
          MALLOCHOST(&(data->_h_px), *(info->_num_atoms) * sizeof(rbmd::Real)));
      CHECK_RUNTIME(
          MALLOCHOST(&(data->_h_py), *(info->_num_atoms) * sizeof(rbmd::Real)));
      CHECK_RUNTIME(
          MALLOCHOST(&(data->_h_pz), *(info->_num_atoms) * sizeof(rbmd::Real)));
      CHECK_RUNTIME(
          MALLOCHOST(&(data->_h_vx), *(info->_num_atoms) * sizeof(rbmd::Real)));
      CHECK_RUNTIME(
          MALLOCHOST(&(data->_h_vy), *(info->_num_atoms) * sizeof(rbmd::Real)));
      CHECK_RUNTIME(
          MALLOCHOST(&(data->_h_vz), *(info->_num_atoms) * sizeof(rbmd::Real)));
      CHECK_RUNTIME(
          MALLOCHOST(&(data->_h_charge), *(info->_num_atoms) * sizeof(rbmd::Real)));
      CHECK_RUNTIME(
          MALLOCHOST(&(data->_h_molecules_id), *(info->_num_atoms) * sizeof(rbmd::Id)));

      CHECK_RUNTIME(
          MALLOCHOST(&(data->_h_flagX), *(info->_num_atoms) * sizeof(rbmd::Id)));
      CHECK_RUNTIME(
          MALLOCHOST(&(data->_h_flagY), *(info->_num_atoms) * sizeof(rbmd::Id)));
      CHECK_RUNTIME(
        MALLOCHOST(&(data->_h_flagZ), *(info->_num_atoms) * sizeof(rbmd::Id)));

  }
  else {
    Logger::Instance().error("\033[31m Unsupported atom_style: {}\033[0m", atom_style );
    exit(EXIT_FAILURE); //
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
          //
          if (line.empty() || line[0] == '#') {
            _line_start = &_mapped_memory[_locate + 1];
            continue;
          }

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
            // image flags
            if (!(iss >> _md_data._structure_data->_h_flagX[index]))
              _md_data._structure_data->_h_flagX[index] = 0;
            if (!(iss >> _md_data._structure_data->_h_flagY[index]))
              _md_data._structure_data->_h_flagY[index] = 0;
            if (!(iss >> _md_data._structure_data->_h_flagZ[index]))
              _md_data._structure_data->_h_flagZ[index] = 0;

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
          //
          if (line.empty() || line[0] == '#') {
            _line_start = &_mapped_memory[_locate + 1];
            continue;
          }

          std::istringstream iss(line);
          if (rbmd::IsLegalLine(line)) {
            iss >> atom_id;
            auto index = atom_id - 1;
            ids[index] = atom_id - 1;
            iss >> atom_type >> data->_h_charge[index];
            iss >> data->_h_px[index] >> data->_h_py[index] >>
                data->_h_pz[index];

            //image flags
            if (!(iss >> data->_h_flagX[index])) data->_h_flagX[index] = 0;
            if (!(iss >> data->_h_flagY[index])) data->_h_flagY[index] = 0;
            if (!(iss >> data->_h_flagZ[index])) data->_h_flagZ[index] = 0;

            types[index] = atom_type - 1;
            ++num;
            // std::cout << atom_id << " " << types[index] << " " <<
            // data->_h_charge[index]  << " " << data->_h_px[index] << " " <<
            // data->_h_py[index] << " " << data->_h_pz[index] << std::endl;
          }
          _line_start = &_mapped_memory[_locate];
        }
      }
    }
    else if ("full" == atom_style)
    {
        rbmd::Id molecules_id;
        auto& full_structure_data = _md_data._structure_data;
        FullStructureData* data =
            dynamic_cast<FullStructureData*>(full_structure_data.get());

        for (auto num = 0; _locate < _file_size && num < atoms_num; ++_locate) {
            if (_mapped_memory[_locate] == '\n') {
                auto line = std::string(_line_start, &_mapped_memory[_locate]);
                //
                if (line.empty() || line[0] == '#') {
                  _line_start = &_mapped_memory[_locate + 1];
                  continue;
                }

                std::istringstream iss(line);
                if (rbmd::IsLegalLine(line))
                {
                    iss >> atom_id;
                    auto index = atom_id - 1;
                    ids[index] = atom_id - 1;

                    iss >> molecules_id >> atom_type >> data->_h_charge[index];
                    iss >> data->_h_px[index] >> data->_h_py[index] >> data->_h_pz[index];

                  // image flags
                  if (!(iss >> data->_h_flagX[index])) data->_h_flagX[index] = 0;
                  if (!(iss >> data->_h_flagY[index])) data->_h_flagY[index] = 0;
                  if (!(iss >> data->_h_flagZ[index])) data->_h_flagZ[index] = 0;

                    types[index] = atom_type - 1;
                    data->_h_molecules_id[index] = molecules_id - 1;
                    MolecularMapInsert(data->_h_molecules_id[index], ids[index]);
                    AtomsMapInsert(types[index], ids[index]);
                    AtomstoMolecular(ids[index], data->_h_molecules_id[index]);
                    ++num;
                    // std::cout << atom_id << " " << data->_h_molecules_id[index] << " " << types[index] << " " <<
                    // data->_h_charge[index]  << " " << data->_h_px[index] << " " <<
                    // data->_h_py[index] << " " << data->_h_pz[index] << std::endl;
                }
                _line_start = &_mapped_memory[_locate];
            }
        }
      SetMolecularGroup();
    }

  } catch (const std::exception& e) {
    // log
    return -1;
  }

  return 0;
}

void AtomicReader::MolecularMapInsert(const rbmd::Id& key, const rbmd::Id& value)
{
    auto it = _molecular_map.find(key);
    if (it != _molecular_map.end())
    {
        it->second.push_back(value);
    }
    else
    {
        _molecular_map.insert(std::make_pair(key, std::vector<rbmd::Id>{ value }));
    }
}

void AtomicReader::AtomsMapInsert(const rbmd::Id& key, const rbmd::Id& value)
{
    auto it = _atoms_map.find(key);
    if (it != _atoms_map.end())
    {
        it->second.push_back(value);
    }
    else
    {
        _atoms_map.insert(std::make_pair(key, std::vector<rbmd::Id>{ value }));
    }

}

void AtomicReader::AtomstoMolecular(const rbmd::Id& key, const rbmd::Id& value)
{
    auto it = _atom_to_molecular_map.find(key);
    if (it != _atom_to_molecular_map.end())
    {
        it->second = value;
    }
    else
    {
        _atom_to_molecular_map.insert(std::make_pair(key, value));
    }
}

void AtomicReader::SetMolecularGroup()
{
    //
    std::vector<std::vector<rbmd::Id>> atoms_gro;
    auto& ids = _md_data._structure_data->_h_atoms_id;
    auto& info = _md_data._structure_info_data;
    for (int index=0;index< *(info->_num_atoms);index++)
    {
        const auto atom_id = ids[index];
        auto molecular_id = _atom_to_molecular_map[atom_id];
        auto atoms_vec = _molecular_map[molecular_id];
        atoms_gro.emplace_back(atoms_vec);
    }

    std::vector<rbmd::Id> atoms_vec_gro;
    std::vector<rbmd::Id> countVector;
    for (const std::vector<rbmd::Id>& innerVector : atoms_gro)
    {
        //flattening
        atoms_vec_gro.insert(atoms_vec_gro.end(), innerVector.begin(), innerVector.end());
        //countVector
        countVector.push_back(innerVector.size());
    }
    auto& full_structure_data = _md_data._structure_data;
    FullStructureData* data =
        dynamic_cast<FullStructureData*>(full_structure_data.get());

    CHECK_RUNTIME(
        MALLOCHOST(&(data->_h_atoms_vec_gro), atoms_vec_gro.size() * sizeof(rbmd::Id)));
    CHECK_RUNTIME(
        MALLOCHOST(&(data->_h_count_vector), countVector.size() * sizeof(rbmd::Id)));
    CHECK_RUNTIME(
     MALLOCHOST(&(data->_h_atoms_offset), (countVector.size()+1)* sizeof(rbmd::Id)));

    memcpy(data->_h_atoms_vec_gro, atoms_vec_gro.data(), atoms_vec_gro.size() * sizeof(rbmd::Id));
    memcpy(data->_h_count_vector, countVector.data(), countVector.size() * sizeof(rbmd::Id));

    std::vector<rbmd::Id> cumulative_offsets;
    cumulative_offsets.push_back(0);
    for (int i = 0; i < countVector.size(); ++i)
    {
      cumulative_offsets.push_back(cumulative_offsets.back() + countVector[i]); //cpu
    }
    memcpy(data->_h_atoms_offset, cumulative_offsets.data(), cumulative_offsets.size() * sizeof(rbmd::Id));

    data->_num_atoms_vec_gro = atoms_vec_gro.size();
    data->_num_count_vector = countVector.size();
    data->_num_atoms_offset = cumulative_offsets.size();
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

int AtomicReader::ReadBond(const rbmd::Id& num_bonds) {
    try {
        auto& full_structure_data = _md_data._structure_data;
        FullStructureData* data = dynamic_cast<FullStructureData*>(full_structure_data.get());
        auto& bond_type = data->_h_bond_type;
        auto& bond_id0 = data->_h_bond_id0;
        auto& bond_id1 = data->_h_bond_id1;
        CHECK_RUNTIME(MALLOCHOST(&bond_type, num_bonds * sizeof(rbmd::Id)));
        CHECK_RUNTIME(MALLOCHOST(&bond_id0, num_bonds * sizeof(rbmd::Id)));
        CHECK_RUNTIME(MALLOCHOST(&bond_id1, num_bonds * sizeof(rbmd::Id)));
        rbmd::Id bound_id_value;
        rbmd::Id bound_type_value;
        rbmd::Id bond_id0_value;
        rbmd::Id bond_id1_value;

        _line_start = &_mapped_memory[_locate];
        for (auto num = 0; _locate < _file_size && num < num_bonds; ++_locate)
        {
            if (_mapped_memory[_locate] == '\n')
            {
              auto line = std::string(_line_start, &_mapped_memory[_locate]);
              // 跳过空行和注释行
              if (line.empty() || line[0] == '#') {
                _line_start = &_mapped_memory[_locate + 1];
                continue;
              }

              std::istringstream iss(line);
              if(rbmd::IsLegalLine(line))
                {
                    iss >> bound_id_value >> bound_type_value >> bond_id0_value >> bond_id1_value;
                    //std::cout << bound_id_value << " "<<  bound_type_value << " " << bond_id0_value << " " << bond_id1_value << std::endl;
                    bond_type[bound_id_value - 1] = bound_type_value - 1;
                    bond_id0[bound_id_value - 1] = bond_id0_value -1;
                    bond_id1[bound_id_value - 1] = bond_id1_value -1;
                    ++num;

                    auto bond_id1 = bond_id0_value - 1;
                    auto bond_id2 = bond_id1_value - 1;
                    //special
                    _special_map.insert(std::make_pair(bond_id1, bond_id2));
                    _special_map.insert(std::make_pair(bond_id2, bond_id1));

                      // record 1-2 connections (有序对，避免重复)
                      auto pair = ordered_pair(bond_id1, bond_id2);
                      data->special_pairs_12.push_back(pair);
                }
                _line_start = &_mapped_memory[_locate];
            }
        }
        SetSpecialBonds();
    }
    catch (const std::exception& e) {
        // log
        return -1;
    }

    return 0;
}

int AtomicReader::ReadAngle(const rbmd::Id& num_angles)
{
    try {
        auto& full_structure_data = _md_data._structure_data;
        FullStructureData* data = dynamic_cast<FullStructureData*>(full_structure_data.get());
        auto& angle_type = data->_h_angle_type;
        auto& angle_id0 = data->_h_angle_id0;
        auto& angle_id1 = data->_h_angle_id1;
        auto& angle_id2 = data->_h_angle_id2;
        auto& angle_id_vec = data->_h_angle_id_vec;
        CHECK_RUNTIME(MALLOCHOST(&angle_type, num_angles * sizeof(rbmd::Id)));
        CHECK_RUNTIME(MALLOCHOST(&angle_id0, num_angles * sizeof(rbmd::Id)));
        CHECK_RUNTIME(MALLOCHOST(&angle_id1, num_angles * sizeof(rbmd::Id)));
        CHECK_RUNTIME(MALLOCHOST(&angle_id2, num_angles * sizeof(rbmd::Id)));
        CHECK_RUNTIME(MALLOCHOST(&angle_id_vec, num_angles * sizeof(Id3))); //TODO:qw:Do we need to "3*" ?
        rbmd::Id angle_id_value;
        rbmd::Id angle_type_value;
        rbmd::Id angle_id0_value;
        rbmd::Id angle_id1_value;
        rbmd::Id angle_id2_value;

        _line_start = &_mapped_memory[_locate];
        for (auto num = 0; _locate < _file_size && num < num_angles; ++_locate)
        {
            if (_mapped_memory[_locate] == '\n')
            {

                auto line = std::string(_line_start, &_mapped_memory[_locate]);
                // 跳过空行和注释行
                if (line.empty() || line[0] == '#') {
                  _line_start = &_mapped_memory[_locate + 1];
                  continue;
                }
                std::istringstream iss(line);

                if (rbmd::IsLegalLine(line))
                {
                    iss >> angle_id_value >>angle_type_value >> angle_id0_value >> angle_id1_value >> angle_id2_value;
                    angle_type[angle_id_value - 1] = angle_type_value - 1;
                    angle_id0[angle_id_value - 1] = angle_id0_value - 1;
                    angle_id1[angle_id_value - 1] = angle_id1_value - 1;
                    angle_id2[angle_id_value - 1] = angle_id2_value - 1;
                    angle_id_vec[angle_id_value - 1].x = angle_id0_value - 1;
                    angle_id_vec[angle_id_value - 1].y = angle_id1_value - 1;
                    angle_id_vec[angle_id_value - 1].z = angle_id2_value - 1;

                    ++num;
                  // record 1-3 connections (first and last atoms)
                  auto pair = ordered_pair(angle_id0_value - 1, angle_id2_value - 1);
                  data->special_pairs_13.push_back(pair);
                    //std::cout << angle_type_value << " " << angle_id0_value << " " << angle_id1_value  << " "  << angle_id2_value << std::endl;

                }
                _line_start = &_mapped_memory[_locate];
            }
        }
    }
    catch (const std::exception& e) {
        // log
        return -1;
    }

    return 0;
}

int AtomicReader::ReadDihedrals(const rbmd::Id& num_dihedrals)
{
    try {
        auto& full_structure_data = _md_data._structure_data;
        FullStructureData* data = dynamic_cast<FullStructureData*>(full_structure_data.get());
        auto& dihedral_type = data->_h_dihedral_type;
        auto& dihedral_id0 = data->_h_dihedral_id0;
        auto& dihedral_id1 = data->_h_dihedral_id1;
        auto& dihedral_id2 = data->_h_dihedral_id2;
        auto& dihedral_id3 = data->_h_dihedral_id3;
        CHECK_RUNTIME(MALLOCHOST(&dihedral_type, num_dihedrals * sizeof(rbmd::Id)));
        CHECK_RUNTIME(MALLOCHOST(&dihedral_id0, num_dihedrals * sizeof(rbmd::Id)));
        CHECK_RUNTIME(MALLOCHOST(&dihedral_id1, num_dihedrals * sizeof(rbmd::Id)));
        CHECK_RUNTIME(MALLOCHOST(&dihedral_id2, num_dihedrals * sizeof(rbmd::Id)));
        CHECK_RUNTIME(MALLOCHOST(&dihedral_id3, num_dihedrals * sizeof(rbmd::Id)));
        rbmd::Id dihedral_id_value;
        rbmd::Id dihedral_type_value;
        rbmd::Id dihedral_id0_value;
        rbmd::Id dihedral_id1_value;
        rbmd::Id dihedral_id2_value;
        rbmd::Id dihedral_id3_value;

        _line_start = &_mapped_memory[_locate];
        for (auto num = 0; _locate < _file_size && num < num_dihedrals; ++_locate)
        {
            if (_mapped_memory[_locate] == '\n')
            {

                auto line = std::string(_line_start, &_mapped_memory[_locate]);
                // 跳过空行和注释行
                if (line.empty() || line[0] == '#') {
                  _line_start = &_mapped_memory[_locate + 1];
                  continue;
                }

                std::istringstream iss(line);
                if (rbmd::IsLegalLine(line))
                {
                    iss >> dihedral_id_value >> dihedral_type_value >> dihedral_id0_value >> dihedral_id1_value >> dihedral_id2_value >> dihedral_id3_value;
                    dihedral_type[dihedral_id_value - 1] = dihedral_type_value - 1;
                    dihedral_id0[dihedral_id_value - 1] = dihedral_id0_value - 1;
                    dihedral_id1[dihedral_id_value - 1] = dihedral_id1_value - 1;
                    dihedral_id2[dihedral_id_value - 1] = dihedral_id2_value - 1;
                    dihedral_id3[dihedral_id_value - 1] = dihedral_id3_value - 1;
                    ++num;

                  // record 1-4 connections (first and last atoms)
                  auto pair = ordered_pair(dihedral_id0_value - 1, dihedral_id3_value - 1);
                  data->special_pairs_14.push_back(pair);
                   //std::cout << dihedral_type_value << " " <<dihedral_id0_value << " " << dihedral_id1_value << " " << dihedral_id2_value << " " << dihedral_id3_value<< std::endl;
                }
                _line_start = &_mapped_memory[_locate];
            }
        }
    }
    catch (const std::exception& e) {
        // log
        return -1;
    }

    return 0;
}

int AtomicReader::ReadImpropers(const rbmd::Id& num_impropers)
{
    try {
        auto& full_structure_data = _md_data._structure_data;
        FullStructureData* data = dynamic_cast<FullStructureData*>(full_structure_data.get());
        auto& improper_type = data->_h_improper_type;
        auto& improper_id0 = data->_h_improper_id0;
        auto& improper_id1 = data->_h_improper_id1;
        auto& improper_id2 = data->_h_improper_id2;
        auto& improper_id3 = data->_h_improper_id3;
        CHECK_RUNTIME(MALLOCHOST(&improper_type, num_impropers * sizeof(rbmd::Id)));
        CHECK_RUNTIME(MALLOCHOST(&improper_id0, num_impropers * sizeof(rbmd::Id)));
        CHECK_RUNTIME(MALLOCHOST(&improper_id1, num_impropers * sizeof(rbmd::Id)));
        CHECK_RUNTIME(MALLOCHOST(&improper_id2, num_impropers * sizeof(rbmd::Id)));
        CHECK_RUNTIME(MALLOCHOST(&improper_id3, num_impropers * sizeof(rbmd::Id)));
        rbmd::Id improper_id_value;
        rbmd::Id improper_type_value;
        rbmd::Id improper_id0_value;
        rbmd::Id improper_id1_value;
        rbmd::Id improper_id2_value;
        rbmd::Id improper_id3_value;

        _line_start = &_mapped_memory[_locate];
        for (auto num = 0; _locate < _file_size && num < num_impropers; ++_locate)
        {
            if (_mapped_memory[_locate] == '\n')
            {

                auto line = std::string(_line_start, &_mapped_memory[_locate]);
                // 跳过空行和注释行
                if (line.empty() || line[0] == '#') {
                  _line_start = &_mapped_memory[_locate + 1];
                  continue;
                }
                std::istringstream iss(line);
                if (rbmd::IsLegalLine(line))
                {
                    iss >> improper_id_value >> improper_type_value >> improper_id0_value >> improper_id1_value >> improper_id2_value >> improper_id3_value;
                    improper_type[improper_id_value - 1] = improper_type_value - 1;
                    improper_id0[improper_id_value - 1] = improper_id0_value - 1;
                    improper_id1[improper_id_value - 1] = improper_id1_value - 1;
                    improper_id2[improper_id_value - 1] = improper_id2_value - 1;
                    improper_id3[improper_id_value - 1] = improper_id3_value - 1;
                    ++num;

                   //std::cout << improper_type_value << " " <<improper_id0_value << " " << improper_id1_value << " " << improper_id2_value << " " << improper_id3_value<< std::endl;
                }
                _line_start = &_mapped_memory[_locate];
            }
        }
    }
    catch (const std::exception& e) {
        // log
        return -1;
    }

    return 0;
}

void AtomicReader::SetSpecialBonds0()
{
    auto special_bonds = DataManager::getInstance().getConfigData()->
  GetArray<rbmd::Real>("special_bonds", "hyper_parameters", "extend");

    auto& full_structure_data = _md_data._structure_data;
    FullStructureData* data = dynamic_cast<FullStructureData*>(full_structure_data.get());
    auto& weights = data->_h_special_weights;
    auto& ids = data->_h_special_ids;
    auto& offsets = data->_h_special_offsets;
    auto& special_offset_count = data->_h_special_offset_count;

    std::vector<rbmd::Real> special_weights;
    std::vector<rbmd::Id> special_ids;
    std::vector<rbmd::Id> special_offsets;
    auto& ids_atoms = _md_data._structure_data->_h_atoms_id;
    auto& info = _md_data._structure_info_data;

    for (int i =0;i< *(info->_num_atoms);i++)
    {
        auto atoms_id = ids_atoms[i];
        //non-bond
        if (_special_map.find(atoms_id) == _special_map.end())
        {
            special_weights.push_back(1.0);
            special_ids.push_back(atoms_id);
            special_offsets.push_back(1);
            continue;
        }

        //bond
        rbmd::Id offset = 0;
        auto link_0 = _special_map.equal_range(atoms_id);
        for (auto it0 = link_0.first; it0 != link_0.second; ++it0)
        {
            //1-2 weight
            int key_1 = it0->second;
            special_weights.push_back(special_bonds[0]); // 1-2 weight
            special_ids.push_back(key_1);
            offset++;

            if (_special_map.find(key_1) == _special_map.end())
                continue;

            // 1-3 weight
            auto link_1 = _special_map.equal_range(key_1);
            for (auto it1 = link_1.first; it1 != link_1.second; ++it1)
            {
                auto key_2 = it1->second;
                if (atoms_id == key_2)
                    continue;

                rbmd::Real weight_1_3 = special_bonds[1];
                special_weights.push_back(weight_1_3); // 1-3 weight
                special_ids.push_back(key_2);
                offset++;

                if (_special_map.find(key_2) == _special_map.end())
                    continue;

                //1-4 weight
                auto link_2 = _special_map.equal_range(key_2);
                for (auto it2 = link_2.first; it2 != link_2.second; ++it2)
                {
                    auto key_3 = it2->second;
                    if (key_1 == key_3)
                        continue;

                    rbmd::Real weight_1_4 = special_bonds[2];
                    special_weights.push_back(weight_1_4 ); // 1-4 weight
                    special_ids.push_back(key_3);
                    offset++;
                }
            }
        }

        special_offsets.push_back(offset);
    }

    CHECK_RUNTIME(MALLOCHOST(&weights, special_weights.size() * sizeof(rbmd::Real)));
    CHECK_RUNTIME(MALLOCHOST(&ids, special_ids.size() * sizeof(rbmd::Id)));
    CHECK_RUNTIME(MALLOCHOST(&offsets, (special_offsets.size()+1) * sizeof(rbmd::Id)));
    CHECK_RUNTIME(MALLOCHOST(&special_offset_count, special_offsets.size() * sizeof(rbmd::Id)));

    memcpy(weights, special_weights.data(), special_weights.size() * sizeof(rbmd::Real));
    memcpy(ids, special_ids.data(), special_ids.size() * sizeof(rbmd::Id));
    memcpy(special_offset_count, special_offsets.data(), special_offsets.size() * sizeof(rbmd::Id));

    //cpu:  sum  prefix
    std::vector<rbmd::Id> cumulative_offsets;
    cumulative_offsets.push_back(0); //

    for (size_t i = 0; i < special_offsets.size(); ++i)
    {
      cumulative_offsets.push_back(cumulative_offsets.back() + special_offsets[i]);
    }
    memcpy(offsets, cumulative_offsets.data(), cumulative_offsets.size() * sizeof(rbmd::Id));

    data->_num_special_weights = special_weights.size();
    data->_num_special_ids = special_ids.size();
    data->_num_special_offset_count = special_offsets.size();
    data->_num_special_offsets =cumulative_offsets.size() ;
}

void AtomicReader::SetSpecialBonds()
{
    // 1. 获取参数
    auto special_bonds = DataManager::getInstance().getConfigData()->
        GetArray<rbmd::Real>("special_bonds", "hyper_parameters", "extend");
    // special_bonds[0]: 1-2 weight, [1]: 1-3 weight, [2]: 1-4 weight

    auto& full_structure_data = _md_data._structure_data;
    FullStructureData* data = dynamic_cast<FullStructureData*>(full_structure_data.get());
    auto& weights = data->_h_special_weights;
    auto& ids = data->_h_special_ids;
    auto& offsets = data->_h_special_offsets;
    auto& special_offset_count = data->_h_special_offset_count;

    // 2. 临时容器
    std::vector<rbmd::Real> special_weights_vec;
    std::vector<rbmd::Id> special_ids_vec;
    std::vector<rbmd::Id> special_offsets_vec;

    auto& ids_atoms = data->_h_atoms_id;
    auto num_atoms = *(_md_data._structure_info_data->_num_atoms);

    // 3. 遍历所有原子，寻找拓扑邻居
    for (int i = 0; i < num_atoms; i++)
    {
        auto atom_id = ids_atoms[i];

        // 核心：使用 Map 记录最短拓扑距离，防止重复和权重覆盖
        // key: neighbor_id, value: distance (1, 2, 3)
        std::map<rbmd::Id, int> neighbors_dist;

        // --- Layer 1 (1-2 Bonds) ---
        std::vector<rbmd::Id> layer1;
        auto range1 = _special_map.equal_range(atom_id);
        for (auto it = range1.first; it != range1.second; ++it) {
            rbmd::Id neighbor = it->second;
            // 只要是直接连接，就是 1-2，距离为 1
            if (neighbors_dist.find(neighbor) == neighbors_dist.end()) {
                neighbors_dist[neighbor] = 1;
                layer1.push_back(neighbor);
            }
        }

        // --- Layer 2 (1-3 Angles) ---
        std::vector<rbmd::Id> layer2;
        for (auto n1 : layer1) {
            auto range2 = _special_map.equal_range(n1);
            for (auto it = range2.first; it != range2.second; ++it) {
                rbmd::Id neighbor = it->second;
                if (neighbor == atom_id) continue; // 排除自己

                // 只有之前没出现过的才标记为 1-3 (距离 2)
                if (neighbors_dist.find(neighbor) == neighbors_dist.end()) {
                    neighbors_dist[neighbor] = 2;
                    layer2.push_back(neighbor);
                }
            }
        }

        // --- Layer 3 (1-4 Dihedrals) ---
        for (auto n2 : layer2) {
            auto range3 = _special_map.equal_range(n2);
            for (auto it = range3.first; it != range3.second; ++it) {
                rbmd::Id neighbor = it->second;
                if (neighbor == atom_id) continue; // 排除自己

                // 只有之前没出现过的才标记为 1-4 (距离 3)
                if (neighbors_dist.find(neighbor) == neighbors_dist.end()) {
                    neighbors_dist[neighbor] = 3;
                }
            }
        }

        // 4. 将结果存入 flat vectors
        rbmd::Id count = 0;
        if (neighbors_dist.empty()) {
             // 孤立原子占位，保持格式一致
             special_weights_vec.push_back(1.0);
             special_ids_vec.push_back(atom_id);
             count = 1;
        } else {
            for (auto const& [neigh_id, dist] : neighbors_dist) {
                special_ids_vec.push_back(neigh_id);
                // dist 是 1, 2, 3，对应 index 0, 1, 2
                special_weights_vec.push_back(special_bonds[dist - 1]);
                count++;
            }
        }
        special_offsets_vec.push_back(count);
    }

    // 5. 分配内存并拷贝
    CHECK_RUNTIME(MALLOCHOST(&weights, special_weights_vec.size() * sizeof(rbmd::Real)));
    CHECK_RUNTIME(MALLOCHOST(&ids, special_ids_vec.size() * sizeof(rbmd::Id)));
    CHECK_RUNTIME(MALLOCHOST(&offsets, (special_offsets_vec.size()+1) * sizeof(rbmd::Id)));
    CHECK_RUNTIME(MALLOCHOST(&special_offset_count, special_offsets_vec.size() * sizeof(rbmd::Id)));

    memcpy(weights, special_weights_vec.data(), special_weights_vec.size() * sizeof(rbmd::Real));
    memcpy(ids, special_ids_vec.data(), special_ids_vec.size() * sizeof(rbmd::Id));
    memcpy(special_offset_count, special_offsets_vec.data(), special_offsets_vec.size() * sizeof(rbmd::Id));

    // 6. 计算前缀和 (offsets)
    std::vector<rbmd::Id> cumulative_offsets;
    cumulative_offsets.push_back(0);
    for (size_t i = 0; i < special_offsets_vec.size(); ++i)
    {
      cumulative_offsets.push_back(cumulative_offsets.back() + special_offsets_vec[i]);
    }
    memcpy(offsets, cumulative_offsets.data(), cumulative_offsets.size() * sizeof(rbmd::Id));

    // 7. 更新计数
    data->_num_special_weights = special_weights_vec.size();
    data->_num_special_ids = special_ids_vec.size();
    data->_num_special_offset_count = special_offsets_vec.size();
    data->_num_special_offsets = cumulative_offsets.size();

    std::cout << " SetSpecialBonds (BFS) completed." << std::endl;
}