#pragma once
#include "common/types.h"
#include "neighbor_list/include/linked_cell.h"

namespace op {
template <typename DEVICE>
struct InitlizeCellOp {
  void operator()(LinkedCellDeviceDataPtr* linked_cell, Box* box, Cell* cells,
                  rbmd::Id total_cells);
};

template <typename DEVICE>
struct AssignAtomsToCellOp {
  void operator()(rbmd::Real* px, rbmd::Real* py, rbmd::Real* pz, Box* d_box,
                  LinkedCellDeviceDataPtr* linked_cell, Cell* cells,
                  rbmd::Id* per_atom_cell_id, rbmd::Id total_atoms_num);
};

template <typename DEVICE>
struct ComputeCellRangesIndicesOp {
  void operator()(rbmd::Id* sorted_cell_index,
                  rbmd::Id* d_in_atom_list_start_index,
                  rbmd::Id* d_in_atom_list_end_index, rbmd::Id num_atoms);
};
} // namespace op
