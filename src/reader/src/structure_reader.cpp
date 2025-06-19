#include "../include/structure_reader.h"

#include "rbmd_define.h"

#include <sstream>

#include "../Utilities/string_util.h"
#include "model/md_data.h"
#include "data_manager.h"
#include "output/include/Logger.hpp"
StructureReder::StructureReder(const std::string& filePath, MDData& data)
    : MmapReader(filePath), _md_data(data) {}

int StructureReder::Execute() {
  if (-1 == MmapReader::Execute()) {
    // log
    return -1;
  }

  if (-1 == ReadHeader()) {
    // log
    return -1;
  }

  if (-1 == ReadForceField()) {
    // log
    return -1;
  }

  AllocateDataSpace();
  if (-1 == ReadData()) {
    // log
    return -1;
  }

  return 0;
}

int StructureReder::ReadHeader() {
  try {
    auto& info = _md_data._structure_info_data;
    CHECK_RUNTIME(MALLOCHOST(&(info->_num_atoms), sizeof(rbmd::Id)));
    CHECK_RUNTIME(MALLOCHOST(&(info->_num_bonds), sizeof(rbmd::Id)));
    CHECK_RUNTIME(MALLOCHOST(&(info->_num_angles), sizeof(rbmd::Id)));
    CHECK_RUNTIME(MALLOCHOST(&(info->_num_dihedrals), sizeof(rbmd::Id)));
    CHECK_RUNTIME(MALLOCHOST(&(info->_num_impropers), sizeof(rbmd::Id)));
    CHECK_RUNTIME(MALLOCHOST(&(info->_num_atoms_type), sizeof(rbmd::Id)));
    CHECK_RUNTIME(MALLOCHOST(&(info->_num_bounds_type), sizeof(rbmd::Id)));
    CHECK_RUNTIME(MALLOCHOST(&(info->_num_angles_type), sizeof(rbmd::Id)));
    CHECK_RUNTIME(MALLOCHOST(&(info->_num_dihedrals_type), sizeof(rbmd::Id)));
    CHECK_RUNTIME(MALLOCHOST(&(info->_num_impropers_type), sizeof(rbmd::Id)));
    CHECK_RUNTIME(MALLOCHOST(&(info->_range), sizeof(rbmd::Range)));
    auto& box = _md_data._box;
    rbmd::Real coord_min[3];
    rbmd::Real coord_max[3];

    for (; _locate < _file_size; ++_locate) {
      if (_mapped_memory[_locate] == '\n') {
        auto line = std::string(_line_start, &_mapped_memory[_locate]);
        std::istringstream iss(line);
        if (rbmd::IsLegalLine(line)) {
          if (line.find("atoms") != std::string::npos) {
            iss >> *(info->_num_atoms);
            Logger::Instance().info("read {} atoms", *(info->_num_atoms));
          } else if (line.find("bonds") != std::string::npos) {
            iss >> *(info->_num_bonds);
            Logger::Instance().info("read {} bonds", *(info->_num_bonds));
          } else if (line.find("angles") != std::string::npos) {
            iss >> *(info->_num_angles);
            Logger::Instance().info("read {} angles", *(info->_num_angles));
          } else if (line.find("dihedrals") != std::string::npos) {
            iss >> *(info->_num_dihedrals);
            Logger::Instance().info("read {} dihedrals",*(info->_num_dihedrals));
          } else if (line.find("impropers") != std::string::npos) {
            iss >> *(info->_num_impropers);
            Logger::Instance().info("read {} impropers",*(info->_num_impropers));
          } else if (line.find("atom types") != std::string::npos) {
            iss >> *(info->_num_atoms_type);
            Logger::Instance().info("read {} atom types",*(info->_num_atoms_type));
          } else if (line.find("bond types") != std::string::npos) {
            iss >> *(info->_num_bounds_type);
            Logger::Instance().info("read {} bond types",*(info->_num_bounds_type));
          } else if (line.find("angle types") != std::string::npos) {
            iss >> *(info->_num_angles_type);
            Logger::Instance().info("read {} angle types",*(info->_num_angles_type));
          }
          else if (line.find("dihedral types") != std::string::npos) {
              iss >> *(info->_num_dihedrals_type);
            Logger::Instance().info("read {} dihedral types",*(info->_num_dihedrals_type));
          } else if (line.find("improper types") != std::string::npos) {
            iss >> *(info->_num_impropers_type);
            Logger::Instance().info("read {} improper types",*(info->_num_impropers_type));
          } else if (line.find("xlo xhi") != std::string::npos) {
            iss >> coord_min[0] >> coord_max[0];
            // std::cout << coord_min[0] << " " << coord_max[0] << "xlo xhi" <<
            // std::endl;
            (*info->_range)[0][0] = coord_min[0];
            (*info->_range)[0][1] = coord_max[0];
          } else if (line.find("ylo yhi") != std::string::npos) {
            iss >> coord_min[1] >> coord_max[1];
            // std::cout << coord_min[1] << " " << coord_max[1] << "ylo yhi" <<
            // std::endl;
            (*info->_range)[1][0] = coord_min[1];
            (*info->_range)[1][1] = coord_max[1];
          } else if (line.find("zlo zhi") != std::string::npos) {
            iss >> coord_min[2] >> coord_max[2];
            // std::cout << coord_min[2] << " " << coord_max[2] << "zlo zhi" <<
            // std::endl;
            (*info->_range)[2][0] = coord_min[2];
            (*info->_range)[2][1] = coord_max[2];
            bool pbc[3] = {1, 1, 1};
            box->Setup(box->_type, coord_min, coord_max, pbc);
            // 计算盒子边长
            double length_x = coord_max[0] - coord_min[0];
            double length_y = coord_max[1] - coord_min[1];
            double length_z = coord_max[2] - coord_min[2];
            std::string box_type = "ORTHOGONAL";
            if (box->_type == Box::BoxType::TRICLINIC) {
              box_type = "TRICLINIC";
            }
            Logger::Instance().info(
            "Initial Box Configuration:\n"
            "         Type:       {}\n"
            "         Min (Å):    ({:.2f}, {:.2f}, {:.2f})\n"
            "         Max (Å):    ({:.2f}, {:.2f}, {:.2f})\n"
            "         Lengths (Å): ({:.2f}, {:.2f}, {:.2f})\n"
            "         PBC:        x={}, y={}, z={}",
            box_type,
            coord_min[0], coord_min[1], coord_min[2],
            coord_max[0], coord_max[1], coord_max[2],
            length_x, length_y, length_z,
            pbc[0] ? "periodic" : "fixed",
            pbc[1] ? "periodic" : "fixed",
            pbc[2] ? "periodic" : "fixed"
            );
            _line_start = &_mapped_memory[_locate];
            break;
          }
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

int StructureReder::ReadForceField() {
  try {
    auto& info = _md_data._structure_info_data;
    for (; _locate < _file_size; ++_locate) {
      if (_mapped_memory[_locate] == '\n') {
        auto line = std::string(_line_start, &_mapped_memory[_locate]);
        std::istringstream iss(line);

        if (rbmd::IsLegalLine(line)) {
          if (line.find("Masses") != std::string::npos) {
            // std::cout << "Masses" << std::endl;
            ReadMass(*(info->_num_atoms_type));
          } else if (line.find("Pair Coeffs") != std::string::npos) {
            // std::cout << "Pair Coeffs" << std::endl;
            ReadPairCoeffs(*(info->_num_atoms_type));
          } else if (line.find("Bond Coeffs") != std::string::npos) {
            ReadBondCoeffs(*(info->_num_bounds_type));
            // std::cout << "Bond Coeffs" << std::endl;
          } else if (line.find("Angle Coeffs") != std::string::npos) {
            ReadAngleCoeffs(*(info->_num_angles_type));
            // std::cout << "Angle Coeffs" << std::endl;
          } else if (line.find("Dihedral Coeffs") != std::string::npos) {
            ReadDihedralsCoeffs(*(info->_num_dihedrals_type));
          } else if (line.find("Improper Coeffs") != std::string::npos) {
            ReadImproperCoeffs(*(info->_num_impropers_type));
          } else if (line.find("group") != std::string::npos) {
            break;
          } else if (line.find("Atoms") != std::string::npos) {
            break;
          }
        }
      }
    }
  } catch (const std::exception& e) {
    // log
    return -1;
  }

  return 0;
}

int StructureReder::ReadMass(const rbmd::Id& numAtomTypes) {
  try {
    auto force_filed =
        std::dynamic_pointer_cast<ForceFieldData>(_md_data._force_field_data);
    auto& mass = force_filed->_h_mass;
    CHECK_RUNTIME(MALLOCHOST(&mass, numAtomTypes * sizeof(rbmd::Id)));
    rbmd::Id atom_type;
    rbmd::Real value;

    _line_start = &_mapped_memory[_locate];
    for (auto num = 0; _locate < _file_size && num < numAtomTypes; ++_locate) {
      if (_mapped_memory[_locate] == '\n') {
        auto line = std::string(_line_start, &_mapped_memory[_locate]);
        std::istringstream iss(line);
        if (rbmd::IsLegalLine(line)) {
          iss >> atom_type >> value;
          mass[atom_type - 1] = value;
          // std::cout << atom_type << " " << force_filed->_h_mass[atom_type -
          // 1] << std::endl;
          ++num;
        }

        _line_start = &_mapped_memory[_locate];
      }
    }
    auto& mass_1 = force_filed->_h_mass;
    /*std::cout << "mass[0]=" << force_filed->_h_mass[0] << ","
        << "mass[1]=" << force_filed->_h_mass[1] << std::endl;*/
  }
  catch (const std::exception& e) {
    // log
    return -1;
  }

  return 0;
}

int StructureReder::ReadPairCoeffs(const rbmd::Id& numAtomTypes) {
    auto force_style = DataManager::getInstance().getConfigData()->Get<std::string>("type", "hyper_parameters", "force_field");
    if ("CVFF" == force_style) {
        try {
            auto force_filed =
                std::dynamic_pointer_cast<CVFFForceFieldData>(_md_data._force_field_data);
            auto& eps = force_filed->_h_eps;
            auto& sigma = force_filed->_h_sigma;
            CHECK_RUNTIME(MALLOCHOST(&eps, numAtomTypes * sizeof(rbmd::Id)));
            CHECK_RUNTIME(MALLOCHOST(&sigma, numAtomTypes * sizeof(rbmd::Id)));
            rbmd::Id atom_type;
            rbmd::Real eps_value;
            rbmd::Real sigma_value;

            _line_start = &_mapped_memory[_locate];
            for (auto num = 0; _locate < _file_size && num < numAtomTypes; ++_locate) {
                if (_mapped_memory[_locate] == '\n') {
                    auto line = std::string(_line_start, &_mapped_memory[_locate]);
                    std::istringstream iss(line);
                    if (rbmd::IsLegalLine(line)) {
                        iss >> atom_type >> eps_value >> sigma_value;
                        eps[atom_type - 1] = eps_value;
                        sigma[atom_type - 1] = sigma_value;
                        // std::cout << atom_type << " " << force_filed->_h_eps[atom_type - 1]
                        // << " " << force_filed->_h_sigma[atom_type - 1] << std::endl;
                        ++num;
                    }
                    _line_start = &_mapped_memory[_locate];
                }
            }
        }
        catch (const std::exception& e) {
            // log
            return -1;
        }
    }
    else {
        try {
            auto force_filed =
                std::dynamic_pointer_cast<LJForceFieldData>(_md_data._force_field_data);
            auto& eps = force_filed->_h_eps;
            auto& sigma = force_filed->_h_sigma;
            CHECK_RUNTIME(MALLOCHOST(&eps, numAtomTypes * sizeof(rbmd::Id)));
            CHECK_RUNTIME(MALLOCHOST(&sigma, numAtomTypes * sizeof(rbmd::Id)));
            rbmd::Id atom_type;
            rbmd::Real eps_value;
            rbmd::Real sigma_value;

            _line_start = &_mapped_memory[_locate];
            for (auto num = 0; _locate < _file_size && num < numAtomTypes; ++_locate) {
                if (_mapped_memory[_locate] == '\n') {
                    auto line = std::string(_line_start, &_mapped_memory[_locate]);
                    std::istringstream iss(line);
                    if (rbmd::IsLegalLine(line)) {
                        iss >> atom_type >> eps_value >> sigma_value;
                        eps[atom_type - 1] = eps_value;
                        sigma[atom_type - 1] = sigma_value;
                        // std::cout << atom_type << " " << force_filed->_h_eps[atom_type - 1]
                        // << " " << force_filed->_h_sigma[atom_type - 1] << std::endl;
                        ++num;
                    }
                    _line_start = &_mapped_memory[_locate];
                }
            }
        }
        catch (const std::exception& e) {
            // log
            return -1;
        }
    }
  

  return 0;
}

int StructureReder::ReadBondCoeffs(const rbmd::Id& numBondTypes) {
  const auto& config = DataManager::getInstance().getConfigData();
  if (!config->PathExists({"hyper_parameters", "force_field", "bond_type"})) {
    Logger::Instance().error( "\033[31m Missing 'bond_type' definition "
             "in 'force_field'.\033[0m");
    exit(EXIT_FAILURE); //
  }
  //
  std::string bond_type ="null";
  bond_type = config->Get<std::string>("bond_type", "hyper_parameters", "force_field");
  if (bond_type == "harmonic") {
    try {
      auto force_filed = std::dynamic_pointer_cast<CVFFForceFieldData>(_md_data._force_field_data);
      auto& bond_coeffs_k = force_filed->_h_bond_coeffs_k;
      auto& bond_coeffs_equilibrium = force_filed->_h_bond_coeffs_equilibrium;
      CHECK_RUNTIME(MALLOCHOST(&bond_coeffs_k, numBondTypes * sizeof(rbmd::Real)));
      CHECK_RUNTIME(MALLOCHOST(&bond_coeffs_equilibrium, numBondTypes * sizeof(rbmd::Real)));
      rbmd::Id bound_type;
      rbmd::Real bond_coeffs_k_value;
      rbmd::Real equilibrium_value;

      _line_start = &_mapped_memory[_locate];
      for (auto num = 0; _locate < _file_size && num < numBondTypes; ++_locate)
      {
        if (_mapped_memory[_locate] == '\n')
        {
          auto line = std::string(_line_start,&_mapped_memory[_locate]);
          std::istringstream iss(line);
          if(rbmd::IsLegalLine(line))
          {
            iss >> bound_type >> bond_coeffs_k_value >> equilibrium_value;
            //std::cout << bound_type << " " <<bond_coeffs_k_value << " " << equilibrium_value << std::endl;
            bond_coeffs_k[bound_type - 1] = bond_coeffs_k_value;
            bond_coeffs_equilibrium[bound_type - 1] =equilibrium_value;
            ++num;
          }
          _line_start = &_mapped_memory[_locate];
        }
      }

    } catch (const std::exception& e) {
      // log
      return -1;
    }
  }
  else {
    Logger::Instance().error("\033[31m Unsupported bond_type: {}\033[0m", bond_type );
    exit(EXIT_FAILURE); //
  }
  return 0;
}

int StructureReder::ReadAngleCoeffs(const rbmd::Id& numAngleTypes)
{
  const auto& config = DataManager::getInstance().getConfigData();
  if (!config->PathExists({"hyper_parameters", "force_field", "angle_type"})) {
    Logger::Instance().error( "\033[31m Missing 'angle_type' definition "
         "in 'force_field'.\033[0m");
    exit(EXIT_FAILURE); //
  }
  //
  std::string angle_type ="null";
  angle_type = config->Get<std::string>("angle_type", "hyper_parameters", "force_field");
  if (angle_type == "harmonic") {
    try {
      auto force_filed = std::dynamic_pointer_cast<CVFFForceFieldData>(_md_data._force_field_data);
      auto& angle_coeffs_k = force_filed->_h_angle_coeffs_k;
      auto& angle_coeffs_equilibrium = force_filed->_h_angle_coeffs_equilibrium;
      CHECK_RUNTIME(MALLOCHOST(&angle_coeffs_k, numAngleTypes * sizeof(rbmd::Real)));
      CHECK_RUNTIME(MALLOCHOST(&angle_coeffs_equilibrium, numAngleTypes * sizeof(rbmd::Real)));
      rbmd::Id angle_type;
      rbmd::Real angle_coeffs_k_value;
      rbmd::Real equilibrium_value;

      _line_start = &_mapped_memory[_locate];
      for (auto num = 0; _locate < _file_size && num < numAngleTypes; ++_locate)
      {
        if (_mapped_memory[_locate] == '\n')
        {
          auto line = std::string(_line_start,&_mapped_memory[_locate]);
          std::istringstream iss(line);
          if(rbmd::IsLegalLine(line))
          {
            iss >> angle_type >> angle_coeffs_k_value >>equilibrium_value;
            //std::cout << angle_type << " " <<angle_coeffs_k_value << " " << equilibrium_value << std::endl;
            angle_coeffs_k[angle_type - 1] = angle_coeffs_k_value;
            angle_coeffs_equilibrium[angle_type - 1] = equilibrium_value;
            ++num;
          }
          _line_start = &_mapped_memory[_locate];
        }
      }
    } catch (const std::exception& e) {
      // log
      return -1;
    }
  }
  else {
    Logger::Instance().error("\033[31m Unsupported angle_type: {}\033[0m", angle_type );
    exit(EXIT_FAILURE); //
  }
  return 0;
}

int StructureReder::ReadDihedralsCoeffs(const rbmd::Id& numDihedralsTypes)
{
  const auto& config = DataManager::getInstance().getConfigData();
  if (!config->PathExists({"hyper_parameters", "force_field", "dihedral_type"})) {
    Logger::Instance().error( "\033[31m Missing 'dihedral_type' definition "
     "in 'force_field'.\033[0m");
    exit(EXIT_FAILURE); //
  }

  std::string dihedral_type ="null";
  dihedral_type = config->Get<std::string>("dihedral_type", "hyper_parameters", "force_field");
  if (dihedral_type == "harmonic") {
      try {
          auto force_filed = std::dynamic_pointer_cast<CVFFForceFieldData>(_md_data._force_field_data);
          auto& dihedral_coeffs_k = force_filed->_h_dihedral_coeffs_k;
          auto& dihedral_coeffs_sign = force_filed->_h_dihedral_coeffs_sign;
          auto& dihedral_coeffs_multiplicity = force_filed->_h_dihedral_coeffs_multiplicity;
          CHECK_RUNTIME(MALLOCHOST(&dihedral_coeffs_k, numDihedralsTypes * sizeof(rbmd::Real)));
          CHECK_RUNTIME(MALLOCHOST(&dihedral_coeffs_sign, numDihedralsTypes * sizeof(rbmd::Real)));
          CHECK_RUNTIME(MALLOCHOST(&dihedral_coeffs_multiplicity, numDihedralsTypes * sizeof(rbmd::Real)));
          rbmd::Id dihedral_type_id;
          rbmd::Real dihedral_coeffs_k_value;
          rbmd::Real dihedral_coeffs_sign_value;
          rbmd::Real dihedral_coeffs_multiplicity_value;

          _line_start = &_mapped_memory[_locate];
          for (auto num = 0; _locate < _file_size && num < numDihedralsTypes; ++_locate)
          {
              if (_mapped_memory[_locate] == '\n')
              {
                  auto line = std::string(_line_start, &_mapped_memory[_locate]); std::istringstream iss(line);
                  if (rbmd::IsLegalLine(line))
                  {
                      iss >> dihedral_type_id >> dihedral_coeffs_k_value >> dihedral_coeffs_sign_value >> dihedral_coeffs_multiplicity_value;
                      dihedral_coeffs_k[dihedral_type_id - 1] = dihedral_coeffs_k_value;
                      dihedral_coeffs_sign[dihedral_type_id - 1] = dihedral_coeffs_sign_value;
                      dihedral_coeffs_multiplicity[dihedral_type_id - 1] = dihedral_coeffs_multiplicity_value;
                      //std::cout << dihedral_type_id << " " << dihedral_coeffs_k_value << " " << dihedral_coeffs_sign_value << " " << dihedral_coeffs_multiplicity_value << std::endl;
                      ++num;
                  }
                  _line_start = &_mapped_memory[_locate];
              }
          }
    }
      catch (const std::exception& e) {
        // log
        return -1;
       }
  }
  else if (dihedral_type == "opls") {
    try {
        auto force_filed = std::dynamic_pointer_cast<CVFFForceFieldData>(_md_data._force_field_data);
        auto& dihedral_coeffs_k1 = force_filed->_h_dihedral_coeffs_k1;
        auto& dihedral_coeffs_k2 = force_filed->_h_dihedral_coeffs_k2;
        auto& dihedral_coeffs_k3 = force_filed->_h_dihedral_coeffs_k3;
        auto& dihedral_coeffs_k4 = force_filed->_h_dihedral_coeffs_k4;
        CHECK_RUNTIME(MALLOCHOST(&dihedral_coeffs_k1, numDihedralsTypes * sizeof(rbmd::Real)));
        CHECK_RUNTIME(MALLOCHOST(&dihedral_coeffs_k2, numDihedralsTypes * sizeof(rbmd::Real)));
        CHECK_RUNTIME(MALLOCHOST(&dihedral_coeffs_k3, numDihedralsTypes * sizeof(rbmd::Real)));
        CHECK_RUNTIME(MALLOCHOST(&dihedral_coeffs_k4, numDihedralsTypes * sizeof(rbmd::Real)));
        rbmd::Id dihedral_type_id;
        rbmd::Real dihedral_coeffs_k1_value;
        rbmd::Real dihedral_coeffs_k2_value;
        rbmd::Real dihedral_coeffs_k3_value;
        rbmd::Real dihedral_coeffs_k4_value;

        _line_start = &_mapped_memory[_locate];
        for (auto num = 0; _locate < _file_size && num < numDihedralsTypes; ++_locate)
        {
            if (_mapped_memory[_locate] == '\n')
            {
                auto line = std::string(_line_start, &_mapped_memory[_locate]); std::istringstream iss(line);
                if (rbmd::IsLegalLine(line))
                {
                    iss >> dihedral_type_id >> dihedral_coeffs_k1_value >> dihedral_coeffs_k2_value >> dihedral_coeffs_k3_value >> dihedral_coeffs_k4_value;
                    dihedral_coeffs_k1[dihedral_type_id - 1] = dihedral_coeffs_k1_value;
                    dihedral_coeffs_k2[dihedral_type_id - 1] = dihedral_coeffs_k2_value;
                    dihedral_coeffs_k3[dihedral_type_id - 1] = dihedral_coeffs_k3_value;
                    dihedral_coeffs_k4[dihedral_type_id - 1] = dihedral_coeffs_k4_value;
                    //std::cout << dihedral_type_id << " " << dihedral_coeffs_k_value << " " << dihedral_coeffs_sign_value << " " << dihedral_coeffs_multiplicity_value << std::endl;
                    ++num;
                }
                _line_start = &_mapped_memory[_locate];
            }
        }
    }
    catch (const std::exception& e) {
        // log
        return -1;
    }
  }
  else {
    Logger::Instance().error("\033[31m Unsupported dihedral_type: {}\033[0m", dihedral_type );
    exit(EXIT_FAILURE); //
  }

  return 0;
}

int StructureReder::ReadImproperCoeffs(const rbmd::Id& numImproperTypes)
{
  const auto& config = DataManager::getInstance().getConfigData();
  if (!config->PathExists({"hyper_parameters", "force_field", "improper_type"})) {
    Logger::Instance().error( "\033[31m Missing 'improper_type' definition "
      "in 'force_field'.\033[0m");
    exit(EXIT_FAILURE); //
  }

  std::string improper_type ="null";
  improper_type = DataManager::getInstance().getConfigData()->Get
    <std::string>("improper_type", "hyper_parameters", "force_field");

  if (improper_type == "harmonic") {
        try {
        auto force_filed = std::dynamic_pointer_cast<CVFFForceFieldData>(_md_data._force_field_data);
        auto& improper_coeffs_k = force_filed->_h_improper_coeffs_k;
        auto& improper_coeffs_degree = force_filed->_h_improper_coeffs_degree;

        CHECK_RUNTIME(MALLOCHOST(&improper_coeffs_k, numImproperTypes * sizeof(rbmd::Real)));
        CHECK_RUNTIME(MALLOCHOST(&improper_coeffs_degree, numImproperTypes * sizeof(rbmd::Real)));

        rbmd::Id improper_type_id;
        rbmd::Real improper_coeffs_k_value;
        rbmd::Real improper_coeffs_degree_value;

        _line_start = &_mapped_memory[_locate];
        for (auto num = 0; _locate < _file_size && num < numImproperTypes; ++_locate)
        {
            if (_mapped_memory[_locate] == '\n')
            {
                auto line = std::string(_line_start, &_mapped_memory[_locate]); std::istringstream iss(line);
                if (rbmd::IsLegalLine(line))
                {
                    iss >> improper_type_id >> improper_coeffs_k_value >> improper_coeffs_degree_value;
                    improper_coeffs_k[improper_type_id - 1] = improper_coeffs_k_value;
                    improper_coeffs_degree[improper_type_id - 1] = improper_coeffs_degree_value;
                    //std::cout << improper_type_id << " " << improper_coeffs_k_value << " " <<improper_coeffs_degree_value << std::endl;
                    ++num;
                }
                _line_start = &_mapped_memory[_locate];
            }
        }
    }
    catch (const std::exception& e) {
        // log
        return -1;
    }
  }
  else if (improper_type == "cvff") {
    try {
        auto force_filed = std::dynamic_pointer_cast<CVFFForceFieldData>(_md_data._force_field_data);
        auto& improper_coeffs_k = force_filed->_h_improper_coeffs_k;
        auto& improper_coeffs_d = force_filed->_h_improper_coeffs_d;
        auto& improper_coeffs_n = force_filed->_h_improper_coeffs_n;
        CHECK_RUNTIME(MALLOCHOST(&improper_coeffs_k, numImproperTypes * sizeof(rbmd::Real)));
        CHECK_RUNTIME(MALLOCHOST(&improper_coeffs_d, numImproperTypes * sizeof(rbmd::Id)));
        CHECK_RUNTIME(MALLOCHOST(&improper_coeffs_n, numImproperTypes * sizeof(rbmd::Id)));
        rbmd::Id improper_type_id;
        rbmd::Real improper_coeffs_k_value;
        rbmd::Id improper_coeffs_d_value;
        rbmd::Id improper_coeffs_n_value;

        _line_start = &_mapped_memory[_locate];
        for (auto num = 0; _locate < _file_size && num < numImproperTypes; ++_locate)
        {
            if (_mapped_memory[_locate] == '\n')
            {
                auto line = std::string(_line_start, &_mapped_memory[_locate]); std::istringstream iss(line);
                if (rbmd::IsLegalLine(line))
                {
                    iss >> improper_type_id >> improper_coeffs_k_value >> improper_coeffs_d_value >> improper_coeffs_n_value;
                    improper_coeffs_k[improper_type_id - 1] = improper_coeffs_k_value;
                    improper_coeffs_d[improper_type_id - 1] = improper_coeffs_d_value;
                    improper_coeffs_n[improper_type_id - 1] = improper_coeffs_n_value;
                    //std::cout << "improper_type :::" <<improper_type << " " << improper_coeffs_k_value << " " << improper_coeffs_d_value << " " << improper_coeffs_n_value << std::endl;
                    ++num;
                }
                _line_start = &_mapped_memory[_locate];
            }
        }
    }
    catch (const std::exception& e) {
        // log
        return -1;
    }
  }
  else {
    Logger::Instance().error("\033[31m Unsupported improper_type: {}\033[0m", improper_type );
    exit(EXIT_FAILURE); //
  }
  return 0;
}