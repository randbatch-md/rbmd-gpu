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
        } else if (line.find("Bonds") != std::string::npos) {
          //std::cout << "Bonds" << std::endl;
          ReadBond(*(_md_data._structure_info_data->_num_bonds));
        } else if (line.find("Angles") != std::string::npos) {
          //std::cout << "Angles" << std::endl;
          ReadAngle(*(_md_data._structure_info_data->_num_angles));
        } else if (line.find("Dihedrals") != std::string::npos) {
          //std::cout << "Dihedrals" << std::endl;
          ReadDihedrals(*(_md_data._structure_info_data->_num_dihedrals));
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
  auto& md_data = DataManager::getInstance().getMDData();
  auto& info = _md_data._structure_info_data;
  if ("atomic" == atom_style) {
    auto& data = _md_data._structure_data;

    HIP_CHECK(MALLOCHOST(&(data->_h_atoms_id),*(info->_num_atoms) * sizeof(rbmd::Id)));
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
  } else if ("full" == atom_style) {
      auto& full_structure_data = _md_data._structure_data;
      FullStructureData* data =
          dynamic_cast<FullStructureData*>(full_structure_data.get());
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
      HIP_CHECK(
          MALLOCHOST(&(data->_h_molecules_id), *(info->_num_atoms) * sizeof(rbmd::Id)));

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
            ids[index] = atom_id - 1;
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
    else if ("full" == atom_style) {
        auto& full_structure_data = _md_data._structure_data;
        FullStructureData* data =
            dynamic_cast<FullStructureData*>(full_structure_data.get());
        for (auto num = 0; _locate < _file_size && num < atoms_num; ++_locate) {
            if (_mapped_memory[_locate] == '\n') {
                auto line = std::string(_line_start, &_mapped_memory[_locate]);
                std::istringstream iss(line);
                if (rbmd::IsLegalLine(line)) {
                    iss >> atom_id;
                    auto index = atom_id - 1;
                    ids[index] = atom_id - 1;
                    iss >> data->_h_molecules_id[index]  >> types[index] >> data->_h_charge[index];
                    iss >> data->_h_px[index] >> data->_h_py[index] >>
                        data->_h_pz[index];
                    ++num;
                    /*std::cout << atom_id << " " << data->_h_molecules_id[index] << " " << types[index] << " " <<
                    data->_h_charge[index]  << " " << data->_h_px[index] << " " <<
                    data->_h_py[index] << " " << data->_h_pz[index] << std::endl;*/
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

int AtomicReader::ReadBond(const rbmd::Id& atoms_num) {
    try {
        auto& full_structure_data = _md_data._structure_data;
        FullStructureData* data = dynamic_cast<FullStructureData*>(full_structure_data.get());
        auto& bond_type = data->_h_bond_type;
        auto& bond_id0 = data->_h_bond_id0;
        auto& bond_id1 = data->_h_bond_id1;
        HIP_CHECK(MALLOCHOST(&bond_type, atoms_num * sizeof(rbmd::Id)))
        HIP_CHECK(MALLOCHOST(&bond_id0, atoms_num * sizeof(rbmd::Id)))
        HIP_CHECK(MALLOCHOST(&bond_id1, atoms_num * sizeof(rbmd::Id)));
        rbmd::Id bound_id_value;
        rbmd::Id bound_type_value;
        rbmd::Real bond_id0_value;
        rbmd::Real bond_id1_value;

        _line_start = &_mapped_memory[_locate];
        for (auto num = 0; _locate < _file_size && num < atoms_num; ++_locate)
        {
            if (_mapped_memory[_locate] == '\n')
            {
                auto line = std::string(_line_start,
                    &_mapped_memory[_locate]); std::istringstream iss(line); if
                    (rbmd::IsLegalLine(line))
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
                }
                _line_start = &_mapped_memory[_locate];
            }
        }
        //SetSpecialBonds();
    }
    catch (const std::exception& e) {
        // log
        return -1;
    }

    return 0;
}

int AtomicReader::ReadAngle(const rbmd::Id& atoms_num)
{
    try {
        auto& full_structure_data = _md_data._structure_data;
        FullStructureData* data = dynamic_cast<FullStructureData*>(full_structure_data.get());
        auto& angle_type = data->_h_angle_type;
        auto& angle_id0 = data->_h_angle_id0;
        auto& angle_id1 = data->_h_angle_id1;
        auto& angle_id2 = data->_h_angle_id2;
        HIP_CHECK(MALLOCHOST(&angle_type, atoms_num * sizeof(rbmd::Id)))
        HIP_CHECK(MALLOCHOST(&angle_id0, atoms_num * sizeof(rbmd::Id)));
        HIP_CHECK(MALLOCHOST(&angle_id1, atoms_num * sizeof(rbmd::Id)))
        HIP_CHECK(MALLOCHOST(&angle_id2, atoms_num * sizeof(rbmd::Id)));
        rbmd::Id angle_id_value;
        rbmd::Id angle_type_value;
        rbmd::Real angle_id0_value;
        rbmd::Real angle_id1_value;
        rbmd::Real angle_id2_value;

        _line_start = &_mapped_memory[_locate];
        for (auto num = 0; _locate < _file_size && num < atoms_num; ++_locate)
        {
            if (_mapped_memory[_locate] == '\n')
            {

                auto line = std::string(_line_start, &_mapped_memory[_locate]); std::istringstream iss(line);
                if (rbmd::IsLegalLine(line))
                {
                    iss >> angle_id_value >>angle_type_value >> angle_id0_value >> angle_id1_value >> angle_id2_value;
                    angle_type[angle_id_value - 1] = angle_type_value - 1;
                    angle_id0[angle_id_value - 1] = angle_id0_value - 1;
                    angle_id1[angle_id_value - 1] = angle_id1_value - 1;
                    angle_id2[angle_id_value - 1] = angle_id2_value - 1;
                    ++num;
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

int AtomicReader::ReadDihedrals(const rbmd::Id& atoms_num)
{
    try {
        auto& full_structure_data = _md_data._structure_data;
        FullStructureData* data = dynamic_cast<FullStructureData*>(full_structure_data.get());
        auto& dihedral_type = data->_h_dihedral_type;
        auto& dihedral_id0 = data->_h_dihedral_id0;
        auto& dihedral_id1 = data->_h_dihedral_id1;
        auto& dihedral_id2 = data->_h_dihedral_id2;
        auto& dihedral_id3 = data->_h_dihedral_id3;
        HIP_CHECK(MALLOCHOST(&dihedral_type, atoms_num * sizeof(rbmd::Id)));
        HIP_CHECK(MALLOCHOST(&dihedral_id0, atoms_num * sizeof(rbmd::Id)));
        HIP_CHECK(MALLOCHOST(&dihedral_id1, atoms_num * sizeof(rbmd::Id)));
        HIP_CHECK(MALLOCHOST(&dihedral_id2, atoms_num * sizeof(rbmd::Id)));
        HIP_CHECK(MALLOCHOST(&dihedral_id3, atoms_num * sizeof(rbmd::Id)));
        rbmd::Id dihedral_id_value;
        rbmd::Id dihedral_type_value;
        rbmd::Id dihedral_id0_value;
        rbmd::Id dihedral_id1_value;
        rbmd::Id dihedral_id2_value;
        rbmd::Id dihedral_id3_value;

        _line_start = &_mapped_memory[_locate];
        for (auto num = 0; _locate < _file_size && num < atoms_num; ++_locate)
        {
            if (_mapped_memory[_locate] == '\n')
            {

                auto line = std::string(_line_start, &_mapped_memory[_locate]); std::istringstream iss(line);
                if (rbmd::IsLegalLine(line))
                {
                    iss >> dihedral_id_value >> dihedral_type_value >> dihedral_id0_value >> dihedral_id1_value >> dihedral_id2_value >> dihedral_id3_value;
                    dihedral_type[dihedral_id_value - 1] = dihedral_type_value - 1;
                    dihedral_id0[dihedral_id_value - 1] = dihedral_id0_value - 1;
                    dihedral_id1[dihedral_id_value - 1] = dihedral_id1_value - 1;
                    dihedral_id2[dihedral_id_value - 1] = dihedral_id2_value - 1;
                    dihedral_id3[dihedral_id_value - 1] = dihedral_id3_value - 1;
                    ++num;
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


void AtomicReader::SetSpecialBonds()
{
    //auto special_bonds = DataManager::getInstance().getConfigData()->Get<std::vector<rbmd::Real>>("special_bonds", "hyper_parameters", "extend");
    //auto& full_structure_data = _md_data._structure_data;
    //FullStructureData* data = dynamic_cast<FullStructureData*>(full_structure_data.get());
    //auto& weights = data._h_special_weights;
    //auto& ids = data._h_special_ids;
    //auto& offsets = data._h_special_offsets;
    //
    //std::vector<rbmd::Real> special_weights;
    //std::vector<rbmd::Id> special_ids;
    //std::vector<rbmd::Id> special_offsets;
    //auto& ids_atoms = _md_data._structure_data->_h_atoms_id;
    //for (int i =0;i< sizeof(ids_atoms);i++)
    //{
    //    auto atoms_id = ids_atoms[i];
    //    //??§Ô???????
    //    if (_special_map.find(atoms_id) == _special_map.end())
    //    {
    //        special_weights.push_back(1.0);
    //        special_ids.push_back(atoms_id);
    //        special_offsets.push_back(1);
    //        continue;
    //    }

    //    //????????
    //    rbmd::Id offset = 0;
    //    auto link_0 = _special_map.equal_range(atoms_id);
    //    for (auto it0 = link_0.first; it0 != link_0.second; ++it0)
    //    {
    //        //???????
    //        int key_1 = it0->second;
    //        special_weights.push_back(special_bonds[0]);
    //        special_ids.push_back(key_1);
    //        offset++;

    //        if (_special_map.find(key_1) == _special_map.end())
    //            continue;

    //        //????????
    //        auto link_1 = _special_map.equal_range(key_1);
    //        for (auto it1 = link_1.first; it1 != link_1.second; ++it1)
    //        {
    //            auto key_2 = it1->second;
    //            if (atoms_id == key_2)
    //                continue;

    //            special_weights.push_back(special_bonds[1]);
    //            special_ids.push_back(key_2);
    //            offset++;

    //            if (_special_map.find(key_2) == _special_map.end())
    //                continue;

    //            //????????
    //            auto link_2 = _special_map.equal_range(key_2);
    //            for (auto it2 = link_2.first; it2 != link_2.second; ++it2)
    //            {
    //                auto key_3 = it2->second;
    //                if (key_1 == key_3)
    //                    continue;

    //                special_weights.push_back(special_bonds[2]);
    //                special_ids.push_back(key_3);
    //                offset++;
    //            }
    //        }
    //    }

    //    special_offsets.push_back(offset);
    //}
    //HIP_CHECK(MALLOCHOST(&weights, special_weights.size() * sizeof(rbmd::Real)));
    //HIP_CHECK(MALLOCHOST(&ids, special_ids.size() * sizeof(rbmd::Id)));
    //HIP_CHECK(MALLOCHOST(&offsets, special_offsets.size() * sizeof(rbmd::Id)));
    //weights = special_weights.data();
    //ids = special_ids.data();
    //offsets = special_offsets.data();

}