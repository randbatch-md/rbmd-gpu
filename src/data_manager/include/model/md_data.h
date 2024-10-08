#pragma once
#include <memory>

#include "box.h"
#include "common/object.h"
#include "force_field/cvff_force_field_data.h"
#include "force_field/eam_force_field_data.h"
#include "force_field/force_field_data.h"
#include "force_field/lj_force_field_data.h"
#include "structure_data/atoms_structure_data.h"
#include "structure_data/charge_structure_data.h"
#include "structure_data/full_structure_data.h"
#include "structure_data/basic_structure_data.h"
#include "structure_data/structure_data.h"
#include "structure_info_data.h"

class MDData : public Object {
 public:
  /**
   * @brief constructor
   */
  MDData() {
    _structure_data = std::make_shared<ChargeStructureData>();
    _structure_info_data = std::make_shared<StructureInfoData>();
    _force_field_data = std::make_shared<LJForceFieldData>();
    _h_box = std::make_shared<Box>();
  }

  /**
   * @brief check data
   * @return true or false
   */
  bool checkMDData() {
    if (false == _structure_data->checkStructure()) {
      // log
      //_console->error("check structure data failed!");
      return false;
    }

    if (false == _force_field_data->checkForceField()) {
      // log
      //_console->error("check force field data failed!");
      return false;
    }

    return true;
  }

  std::shared_ptr<StructureData> _structure_data;
  std::shared_ptr<StructureInfoData> _structure_info_data;
  std::shared_ptr<ForceFieldData> _force_field_data;
  std::shared_ptr<Box> _h_box;
};
