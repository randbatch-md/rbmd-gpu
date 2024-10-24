

#include "rbl_half_neighbor_list_builder.h"

#include <common/types.h>

#include "half_neighbor_list_op.h"
#include "rbl_half_neighbor_list_op.h"
void RblHalfNeighborListBuilder::EstimateNeighborsList() {
  std::cout << "\033[31mresizing neighbor list array...\033[0m" << std::endl;
  rbmd::Id* d_total_max_neighbor_num;
  CHECK_RUNTIME(MALLOC(&d_total_max_neighbor_num, sizeof(rbmd::Id)));
  CHECK_RUNTIME(MEMCPY(d_total_max_neighbor_num,
                       &(_neighbor_list->_h_total_max_neighbor_num),
                       sizeof(rbmd::Id), H2D));
  op::EstimateHalfNeighborListOp<device::DEVICE_GPU>
      estimate_half_neighbor_list_op;
  estimate_half_neighbor_list_op(
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
      _neighbor_cell_num, true);  // rbl always true
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
rbmd::Id RblHalfNeighborListBuilder::GenerateNeighborsList() {
  GetRblParams();
  CHECK_RUNTIME(
      MEMCPY(_d_should_realloc, &(this->should_realloc), sizeof(rbmd::Id), H2D));
  op::GenerateRblHalfNeighborListOp<device::DEVICE_GPU> generate_rbl_half_neighbor_list_op;
  generate_rbl_half_neighbor_list_op(
      thrust::raw_pointer_cast(_linked_cell->_per_atom_cell_id.data()),
      thrust::raw_pointer_cast(_linked_cell->_in_atom_list_start_index.data()),
      thrust::raw_pointer_cast(_linked_cell->_in_atom_list_end_index.data()),
      _trunc_distance_power_2,
      _linked_cell->_cutoff * _linked_cell->_cutoff - EPSILON,
      _linked_cell->_total_atoms_num,
      thrust::raw_pointer_cast(_device_data->_d_px.data()),
      thrust::raw_pointer_cast(_device_data->_d_py.data()),
      thrust::raw_pointer_cast(_device_data->_d_pz.data()),
      thrust::raw_pointer_cast(
          this->_neighbor_list->_d_max_neighbor_num.data()),
      thrust::raw_pointer_cast(this->_neighbor_list->_start_idx.data()),
      thrust::raw_pointer_cast(this->_neighbor_list->_end_idx.data()),
      thrust::raw_pointer_cast(this->_neighbor_list->_d_neighbors.data()),
      _neighbor_sample_num,
      thrust::raw_pointer_cast(this->_neighbor_list->_d_random_neighbor.data()),
      thrust::raw_pointer_cast(
          this->_neighbor_list->_d_random_neighbor_num.data()),
      _d_box, _d_should_realloc,
      thrust::raw_pointer_cast(_linked_cell->_neighbor_cell.data()),
      _neighbor_cell_num, _selection_frequency);
  CHECK_RUNTIME(
      MEMCPY(&(this->should_realloc), _d_should_realloc, sizeof(rbmd::Id), D2H));
  return this->should_realloc;
}