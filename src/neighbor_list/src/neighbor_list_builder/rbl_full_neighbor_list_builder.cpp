
#include "rbl_full_neighbor_list_builder.h"

#include <thrust/detail/raw_pointer_cast.h>

#include "common/device_types.h"
#include "full_neighbor_list_op.h"
#include "rbl_full_neighbor_list_op.h"

RblFullNeighborListBuilder::RblFullNeighborListBuilder() {
  _r_core = DataManager::getInstance().getConfigData()->Get<rbmd::Real>(
      "r_core", "hyper_parameters", "neighbor");
  // _linked_cell->_cell_count_within_cutoff = static_cast<rbmd::Id>(
  //     std::ceil(static_cast<double>(_linked_cell->_cutoff / _r_core)));
  // _linked_cell->Build(DataManager::getInstance().getMDData()->_h_box.get());
  // _linked_cell->SyncHToD();  // TODO 逻辑有点复杂   是否放到locator里面去更好
  // _linked_cell->InitializeCells();
  _neighbor_sample_num =
      DataManager::getInstance().getConfigData()->Get<rbmd::Id>(
          "neighbor_sample_num", "hyper_parameters", "neighbor");
  if (_neighbor_sample_num <= 0) {
    std::cout << "\033[31mError neighbor_sample_num must be large than 0.\033[0m"
                << std::endl;
    exit(0);
  }
  _trunc_distance_power_2 = _r_core * _r_core - EPSILON;  // for estimate
  _neighbor_list->_d_random_neighbor.resize(
      _linked_cell->_total_atoms_num * _neighbor_sample_num);  // cs neighbor
  _neighbor_list->_d_random_neighbor_num.resize(_linked_cell->_total_atoms_num);
  // neighbor num !
  this->_neighbor_cell_num = (2 * _linked_cell->_cell_count_within_cutoff + 1) *
                             (2 * _linked_cell->_cell_count_within_cutoff + 1) *
                             (2 * _linked_cell->_cell_count_within_cutoff + 1);
  if (_linked_cell->_total_cells < this->_neighbor_cell_num) {
    this->_neighbor_cell_num = _linked_cell->_total_cells;
    std::cout << "\033[31mwarning: The current simulation domain is too small "
                 "for PBC to be effective.\033[0m"
              << std::endl;
    this->FullNeighborListBuilder::ComputeNeighborCellsWithoutPBC();
  } else {
    this->FullNeighborListBuilder::ComputeNeighborCells();
  }
}

void RblFullNeighborListBuilder::GetRblParams() {
#pragma region rbl prarms
  _system_rho =
      _linked_cell->_total_atoms_num /
      CalculateVolume(DataManager::getInstance().getMDData()->_h_box.get());
  rbmd::Real coeff_rcs = 1.0 + (0.05 / _system_rho - 0.05);
  rbmd::Id Id_coeff_rcs = std::round(static_cast<double>(coeff_rcs));
  rbmd::Id rs_num = Id_coeff_rcs * _system_rho *
                        std::ceil(4.0 / 3.0 * M_PI * std::pow(_r_core, 3)) +
                    1;
  rbmd::Id rc_num =
      Id_coeff_rcs * _system_rho *
          std::ceil(4.0 / 3.0 * M_PI * std::pow(_linked_cell->_cutoff, 3)) +
      1;
  auto random_rate = static_cast<rbmd::Real>(_neighbor_sample_num) /
                     static_cast<rbmd::Real>(rc_num - rs_num);
  _selection_frequency = std::floor(
      1.0 / random_rate);  // note：The ceil function rarely samples the desired
                           // number of neighbor_sample_num, while the floor
                           // function may to some extent affect performance.
  _neighbor_list->_selection_frequency = this->_selection_frequency;
#pragma endregion
}

rbmd::Id RblFullNeighborListBuilder::GenerateNeighborsList() {
  GetRblParams();
  CHECK_RUNTIME(
      MEMCPY(_d_should_realloc, &(this->should_realloc), sizeof(rbmd::Id), H2D));
  op::GenerateRblFullNeighborListOp<device::DEVICE_GPU>
      generate_full_neighbor_list_op;
  generate_full_neighbor_list_op(
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

