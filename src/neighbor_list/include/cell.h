#pragma once
#include "common/types.h"
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include "common/rbmd_define.h"

struct Cell {
    /// 当前cell的id TODO：考虑时候需要CellManager
    rbmd::Id _cell_id = 0;
    rbmd::Real ALIGN(ALIGN_SIZE(rbmd::Real, 3)) cell_coord_min[3]{};
    rbmd::Real ALIGN(ALIGN_SIZE(rbmd::Real , 3)) cell_coord_max[3]{};
    /// 当前cell内的atom数量
    rbmd::Id _atoms_count = 0;  // TODO
};