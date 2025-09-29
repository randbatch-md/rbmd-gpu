
#include "../neighbor_list_builder/mace_neighbor_list_builder.h"
#include "common/device_types.h"
#include "common/types.h"
#include "data_manager.h"
#include "mace_neighbor_list_op.h"

MACENeighborListBuilder::MACENeighborListBuilder() {
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
    this->MACENeighborListBuilder::ComputeNeighborCellsWithoutPBC();
    //this->MACENeighborListBuilder::ComputeNeighborCells();
  } else {
    this->MACENeighborListBuilder::ComputeNeighborCells();
  }
  _trunc_distance_power_2 = _linked_cell->_cutoff * _linked_cell->_cutoff - EPSILON;
}

std::shared_ptr<NeighborList> MACENeighborListBuilder::Build() {
  _linked_cell->AssignAtomsToCell();
  _linked_cell->SortAtomsByCellKey();
  _linked_cell->ComputeCellRangesIndices();
  if (should_realloc) {
    // 好像就第一入口调用了  todo move to init?
    EstimateNeighborsList();
  }
  // 索引没有问题
  //if (GenerateNeighborsList()==RBMD_TRUE) {
    EstimateNeighborsList();
    GenerateNeighborsList();
  //}
  return _neighbor_list;
}


std::shared_ptr<NeighborList> MACENeighborListBuilder::Build( rbmd::Real custom_cutoff) {
  _trunc_distance_power_2 = custom_cutoff * custom_cutoff;
  _linked_cell->AssignAtomsToCell();
  _linked_cell->SortAtomsByCellKey();
  _linked_cell->ComputeCellRangesIndices();
  if (should_realloc) {
    EstimateNeighborsList();
  }
    EstimateNeighborsList();
    GenerateNeighborsList();
  return _neighbor_list;
}

void MACENeighborListBuilder::ComputeNeighborCells() {
  _linked_cell->_neighbor_cell.resize(
      (_linked_cell->_total_cells * this->_neighbor_cell_num));
  op::ComputeMACENeighborsOp<device::DEVICE_GPU> compute_mace_neighbors_op;
  compute_mace_neighbors_op(
      _linked_cell->GetDataPtr()->_d_per_dimension_cells,
      thrust::raw_pointer_cast(_linked_cell->_neighbor_cell.data()),
      this->_neighbor_cell_num, _linked_cell->_total_cells,
      _linked_cell->_cell_count_within_cutoff);
}
void MACENeighborListBuilder::ComputeNeighborCellsWithoutPBC() {
  _linked_cell->_neighbor_cell.resize(
      (_linked_cell->_total_cells * this->_neighbor_cell_num));
  op::ComputeMACENeighborsWithoutPBCOp<device::DEVICE_GPU>
      compute_mace_neighbors_without_pbc_op;
  compute_mace_neighbors_without_pbc_op(
      thrust::raw_pointer_cast(_linked_cell->_neighbor_cell.data()),
      this->_neighbor_cell_num, _linked_cell->_total_cells);
}

void MACENeighborListBuilder::EstimateNeighborsList() {
  std::cout << "\033[31mresizing neighbor list array...\033[0m" << std::endl;
  rbmd::Id* d_total_max_neighbor_num;
  _trunc_distance_power_2 = _linked_cell->_cutoff * _linked_cell->_cutoff - EPSILON;
  CHECK_RUNTIME(MALLOC(&d_total_max_neighbor_num, sizeof(rbmd::Id)));
  CHECK_RUNTIME(MEMCPY(d_total_max_neighbor_num,
                       &(_neighbor_list->_h_total_max_neighbor_num),
                       sizeof(rbmd::Id), H2D));
  op::EstimateMACENeighborListOp<device::DEVICE_GPU>
      estimate_mace_neighbor_list_op;
  estimate_mace_neighbor_list_op(
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
      *_box,
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

rbmd::Id MACENeighborListBuilder::GenerateNeighborsList() {
  CHECK_RUNTIME(
      MEMCPY(_d_should_realloc, &(this->should_realloc), sizeof(rbmd::Id), H2D));
  op::GenerateMACENeighborListOp<device::DEVICE_GPU> generate_mace_neighbor_list_op;

  //rbmd::Id sum = thrust::reduce(this->_neighbor_list->_d_neighbor_num.begin(), this->_neighbor_list->_d_neighbor_num.end(), rbmd::Id(0), thrust::plus<rbmd::Id>());
  //std::cout<<"边数为"<<sum<<std::endl;
  rbmd::Id edgelong = this->_neighbor_list->_end_idx.back();
 //edgelong = 38400;
  this->_neighbor_list->_d_neighbors_atoms.resize(edgelong);
  this->_neighbor_list->_d_shiftx.resize(edgelong);
  this->_neighbor_list->_d_shifty.resize(edgelong);
  this->_neighbor_list->_d_shiftz.resize(edgelong);
  this->_neighbor_list->_d_unit_shiftx.resize(edgelong);
  this->_neighbor_list->_d_unit_shifty.resize(edgelong);
  this->_neighbor_list->_d_unit_shiftz.resize(edgelong);
  generate_mace_neighbor_list_op(
  thrust::raw_pointer_cast(_linked_cell->_per_atom_cell_id.data()),
 thrust::raw_pointer_cast(_linked_cell->_in_atom_list_start_index.data()),
 thrust::raw_pointer_cast(_linked_cell->_in_atom_list_end_index.data()),
 _trunc_distance_power_2, _linked_cell->_total_atoms_num,
 thrust::raw_pointer_cast(_device_data->_d_px.data()),
 thrust::raw_pointer_cast(_device_data->_d_py.data()),
 thrust::raw_pointer_cast(_device_data->_d_pz.data()),
 thrust::raw_pointer_cast(
     this->_neighbor_list->_d_max_neighbor_num.data()),
 thrust::raw_pointer_cast(this->_neighbor_list->_start_idx.data()),
 thrust::raw_pointer_cast(this->_neighbor_list->_end_idx.data()),
 thrust::raw_pointer_cast(this->_neighbor_list->_d_neighbors.data()),
 *_box,
 thrust::raw_pointer_cast(this->_neighbor_list->_d_neighbors_atoms.data()),
      thrust::raw_pointer_cast(this->_neighbor_list->_d_unit_shiftx.data()),
      thrust::raw_pointer_cast(this->_neighbor_list->_d_unit_shifty.data()),
      thrust::raw_pointer_cast(this->_neighbor_list->_d_unit_shiftz.data()),
      thrust::raw_pointer_cast(this->_neighbor_list->_d_shiftx.data()),
      thrust::raw_pointer_cast(this->_neighbor_list->_d_shifty.data()),
      thrust::raw_pointer_cast(this->_neighbor_list->_d_shiftz.data()),
 _d_should_realloc,
 thrust::raw_pointer_cast(_linked_cell->_neighbor_cell.data()),
 _neighbor_cell_num
  );
  //std::cout<<this->_neighbor_list->_d_neighbor_num.size()<<std::endl;
  /*std::cout<< _neighbor_list->_d_neighbor_num[100]<< " "<< _neighbor_list->_d_neighbor_num[200]<<" "<< this->_neighbor_list->_d_neighbor_num.back() <<std::endl;
  //_d_neighbor_num保持不变，应为bug
  std::cout<< this->_neighbor_list->_start_idx.front() <<" "<< this->_neighbor_list->_start_idx.back() <<" "<< this->_neighbor_list->_end_idx.front()
      <<" "<<this->_neighbor_list->_end_idx.back()<<std::endl;
  std::cout<<"最后一个连接的"<<this->_neighbor_list->_d_neighbors[this->_neighbor_list->_end_idx.back()-1]<<std::endl;*/
  CHECK_RUNTIME(
      MEMCPY(&(this->should_realloc), _d_should_realloc, sizeof(rbmd::Id), D2H));
  return this->should_realloc;
}
