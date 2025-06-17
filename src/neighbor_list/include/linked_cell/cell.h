#pragma once
#include "common/rbmd_define.h"
#include "common/types.h"

struct Cell {
  ///
  rbmd::Id _cell_id = 0;
  rbmd::Real  cell_coord_min[3]{};
  rbmd::Real  cell_coord_max[3]{};
  ///  number of atoms in current  cell
  rbmd::Id _atoms_count = 0;  // TODO
};