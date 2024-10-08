
#include "neighbor_list_builder/full_neighbor_list_builder.h"

#include "common/device_types.h"
#include "common/types.h"
#include "data_manager.h"
#include "full_neighbor_list_op.h"

FullNeighborListBuilder::FullNeighborListBuilder() {
  this->_neighbor_cell_num =
      (2 * _linked_cell->_cell_count_within_cutoff + 1) *
      (2 * _linked_cell->_cell_count_within_cutoff + 1) *
      (2 * _linked_cell->_cell_count_within_cutoff +
       1);  // TODO not test _cell_count_within_cutoff=2 ... just 1
  this->_neighbor_list =
      std::make_shared<NeighborList>(_linked_cell->_total_atoms_num, false);
  this->_device_data = DataManager::getInstance().getDeviceData();
  if (_linked_cell->_total_cells < this->_neighbor_cell_num) {
    this->_neighbor_cell_num = _linked_cell->_total_cells;
    std::cout << "\033[31mwarning: The current simulation domain is too small "
                 "for PBC to be effective.\033[0m"
              << std::endl;
    this->FullNeighborListBuilder::ComputeNeighborCellsWithoutPBC();
  } else {
    this->FullNeighborListBuilder::ComputeNeighborCells();
  }
  _trunc_distance_power_2 = _linked_cell->_cutoff * _linked_cell->_cutoff - EPSILON;
}

std::shared_ptr<NeighborList> FullNeighborListBuilder::Build() {
  _linked_cell->AssignAtomsToCell();
  _linked_cell->SortAtomsByCellKey();
  _linked_cell->ComputeCellRangesIndices();
  if (should_realloc) {
    // 好像就第一入口调用了  todo move to init?
    EstimateNeighborsList();
  }
  // 索引没有问题
  if (GenerateNeighborsList()==RBMD_TRUE) {
    EstimateNeighborsList();
    GenerateNeighborsList();
  }
  return _neighbor_list;
}

void FullNeighborListBuilder::ComputeNeighborCells() {
  _linked_cell->_neighbor_cell.resize(
      (_linked_cell->_total_cells * this->_neighbor_cell_num));
  op::ComputeFullNeighborsOp<device::DEVICE_GPU> compute_full_neighbors_op;
  compute_full_neighbors_op(
      _linked_cell->GetDataPtr()->_d_per_dimension_cells,
      thrust::raw_pointer_cast(_linked_cell->_neighbor_cell.data()),
      this->_neighbor_cell_num, _linked_cell->_total_cells,
      _linked_cell->_cell_count_within_cutoff);
}
void FullNeighborListBuilder::ComputeNeighborCellsWithoutPBC() {
  _linked_cell->_neighbor_cell.resize(
      (_linked_cell->_total_cells * this->_neighbor_cell_num));
  op::ComputeFullNeighborsWithoutPBCOp<device::DEVICE_GPU>
      compute_full_neighbors_without_pbc_op;
  compute_full_neighbors_without_pbc_op(
      thrust::raw_pointer_cast(_linked_cell->_neighbor_cell.data()),
      this->_neighbor_cell_num, _linked_cell->_total_cells);
}

void FullNeighborListBuilder::EstimateNeighborsList() {
  std::cout << "\033[31mresizing neighbor list array...\033[0m" << std::endl;
  rbmd::Id* d_total_max_neighbor_num;
  CHECK_RUNTIME(MALLOC(&d_total_max_neighbor_num, sizeof(rbmd::Id)));
  CHECK_RUNTIME(MEMCPY(d_total_max_neighbor_num,
                       &(_neighbor_list->_h_total_max_neighbor_num),
                       sizeof(rbmd::Id), H2D));
  op::EstimateFullNeighborListOp<device::DEVICE_GPU>
      estimate_full_neighbor_list_op;
  estimate_full_neighbor_list_op(
      thrust::raw_pointer_cast(_linked_cell->_per_atom_cell_id.data()),
      thrust::raw_pointer_cast(_linked_cell->_in_atom_list_start_index.data()),
      thrust::raw_pointer_cast(_linked_cell->_in_atom_list_end_index.data()),
      _trunc_distance_power_2, _linked_cell->_total_atoms_num,
      thrust::raw_pointer_cast(_device_data->_d_px.data()),
      thrust::raw_pointer_cast(_device_data->_d_py.data()),
      thrust::raw_pointer_cast(_device_data->_d_pz.data()),
      thrust::raw_pointer_cast(this->_neighbor_list->_d_neighbor_num.data()),
      thrust::raw_pointer_cast(
          this->_neighbor_list->_d_max_neighbor_num.data()),
      this->_d_box,
      thrust::raw_pointer_cast(_linked_cell->_neighbor_cell.data()),
      _neighbor_cell_num);
  ReductionSum(
      thrust::raw_pointer_cast(_neighbor_list->_d_max_neighbor_num.data()),
      d_total_max_neighbor_num, _linked_cell->_total_atoms_num);
  CHECK_RUNTIME(MEMCPY(&(_neighbor_list->_h_total_max_neighbor_num),
                       d_total_max_neighbor_num, sizeof(rbmd::Id), D2H));
  CHECK_RUNTIME(FREE(d_total_max_neighbor_num));
  _neighbor_list->_d_neighbors.resize(
      _neighbor_list->_h_total_max_neighbor_num);
  InitNeighborListIndices();
  this->should_realloc = false;
}

rbmd::Id FullNeighborListBuilder::GenerateNeighborsList() {
  CHECK_RUNTIME(
      MEMCPY(_d_should_realloc, &(this->should_realloc), sizeof(rbmd::Id), H2D));
  op::GenerateFullNeighborListOp<device::DEVICE_GPU>
      generate_full_neighbor_list_op;
  generate_full_neighbor_list_op(
      thrust::raw_pointer_cast(_linked_cell->_per_atom_cell_id.data()),
      thrust::raw_pointer_cast(_linked_cell->_in_atom_list_start_index.data()),
      thrust::raw_pointer_cast(_linked_cell->_in_atom_list_end_index.data()),
      _trunc_distance_power_2,
      _linked_cell->_total_atoms_num,
      thrust::raw_pointer_cast(_device_data->_d_px.data()),
      thrust::raw_pointer_cast(_device_data->_d_py.data()),
      thrust::raw_pointer_cast(_device_data->_d_pz.data()),
      thrust::raw_pointer_cast(
          this->_neighbor_list->_d_max_neighbor_num.data()),
      thrust::raw_pointer_cast(this->_neighbor_list->_start_idx.data()),
      thrust::raw_pointer_cast(this->_neighbor_list->_end_idx.data()),
      thrust::raw_pointer_cast(this->_neighbor_list->_d_neighbors.data()),
      _d_box, _d_should_realloc,
      thrust::raw_pointer_cast(_linked_cell->_neighbor_cell.data()),
      _neighbor_cell_num);
  CHECK_RUNTIME(
      MEMCPY(&(this->should_realloc), _d_should_realloc, sizeof(rbmd::Id), D2H));
  return this->should_realloc;
}
