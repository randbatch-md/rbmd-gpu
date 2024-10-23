#include "linked_cell/linked_cell.h"

#include <thrust/sort.h>
#include <hipcub/backend/rocprim/device/device_radix_sort.hpp>
#include <hipcub/backend/rocprim/iterator/counting_input_iterator.hpp>

#include "../common/device_types.h"
#include "../common/types.h"
#include "data_manager.h"
#include "linked_cell_op.h"
#include "model/md_data.h"

LinkedCell::LinkedCell() {
  this->_device_data = DataManager::getInstance().getDeviceData();
  this->_structure_info_data =
      DataManager::getInstance().getMDData()->_structure_info_data;
  this->_config_data = DataManager::getInstance().getConfigData();
  this->_per_atom_cell_id.resize(*(_structure_info_data->_num_atoms));
  // TODO  反序列化
  this->_cutoff =  _config_data->Get<rbmd::Real>("cut_off",
                      "hyper_parameters", "neighbor");
  this->_total_atoms_num =
      *(_structure_info_data->_num_atoms);  // do this because nativate num
  this->_atom_id_to_idx.resize(_total_atoms_num);
}

LinkedCell::~LinkedCell() {
  CHECK_RUNTIME(FREE(this->_linked_cell_device_data_ptr));
}

__host__ void LinkedCell::Build(Box* box) {
  rbmd::Id _cells_number = 1;
  auto per_cell_length = _cutoff / _cell_count_within_cutoff;
  for (int dim = 0; dim < 3; dim++) {
    box->_box_width_as_cell_units[dim] = static_cast<rbmd::Id>(
        floor(static_cast<double>(box->_coord_max[dim] - box->_coord_min[dim]) /
              per_cell_length));
    _per_dimension_cells[dim] = box->_box_width_as_cell_units[dim];
    _cells_number *= _per_dimension_cells[dim];

    const rbmd::Real diff = box->_coord_max[dim] - box->_coord_min[dim];
    _cell_length[dim] = diff / box->_box_width_as_cell_units[dim];

    // Calculate start and end indices for each dimension
    box->_length[dim] = box->_coord_max[dim] - box->_coord_min[dim];
    _cell_length_reciprocal[dim] = 1.0 / _cell_length[dim];
  }
  this->_total_cells = _cells_number;
  this->_in_atom_list_start_index.resize(_cells_number);
  this->_in_atom_list_end_index.resize(_cells_number);
  this->_cells.resize(_cells_number);
  // update device box
  CHECK_RUNTIME(MEMCPY(DataManager::getInstance().getDeviceData()->_d_box,
                       DataManager::getInstance().getMDData()->_h_box.get(),
                       sizeof(Box), H2D));
}

LinkedCellDeviceDataPtr* LinkedCell::GetDataPtr() {
  if (nullptr == this->_linked_cell_device_data_ptr) {
    this->AllocDeviceMemory();
    this->SyncHToD();
  }
  return this->_linked_cell_device_data_ptr;
}

void LinkedCell::InitializeCells() {
  op::InitializeCellOp<device::DEVICE_GPU> initialize_cell_op;
  initialize_cell_op(GetDataPtr(), _device_data->_d_box,
                     thrust::raw_pointer_cast(this->_cells.data()),
                     this->_total_cells);
}

void LinkedCell::AssignAtomsToCell() {
  op::AssignAtomsToCellOp<device::DEVICE_GPU> assign_atoms_to_cell_op;
  assign_atoms_to_cell_op(thrust::raw_pointer_cast(_device_data->_d_px.data()),
                          thrust::raw_pointer_cast(_device_data->_d_py.data()),
                          thrust::raw_pointer_cast(_device_data->_d_pz.data()),
                          _device_data->_d_box, GetDataPtr(),
                          thrust::raw_pointer_cast(this->_cells.data()),
                          thrust::raw_pointer_cast(_per_atom_cell_id.data()),
                          this->_total_atoms_num);
}

void LinkedCell::ComputeCellRangesIndices() {
  op::ComputeCellRangesIndicesOp<device::DEVICE_GPU>
      compute_cell_ranges_indices_op;
  compute_cell_ranges_indices_op(
      thrust::raw_pointer_cast(_per_atom_cell_id.data()),
      thrust::raw_pointer_cast(_in_atom_list_start_index.data()),
      thrust::raw_pointer_cast(_in_atom_list_end_index.data()),
      this->_total_atoms_num);
}

void LinkedCell::SyncHToD() {
  if (nullptr == this->_linked_cell_device_data_ptr) {
    this->AllocDeviceMemory();
  }
  CHECK_RUNTIME(MEMCPY(&_linked_cell_device_data_ptr->_d_total_cells,
                       &this->_total_cells, sizeof(rbmd::Id), H2D));

  CHECK_RUNTIME(MEMCPY(_linked_cell_device_data_ptr->_d_per_dimension_cells,
                       this->_per_dimension_cells, ALIGN_SIZE(rbmd::Id, 3),
                       H2D));

  CHECK_RUNTIME(MEMCPY(_linked_cell_device_data_ptr->_d_cell_length,
                       this->_cell_length, ALIGN_SIZE(rbmd::Real, 3), H2D));

  CHECK_RUNTIME(MEMCPY(_linked_cell_device_data_ptr->_d_cell_length_reciprocal,
                       this->_cell_length_reciprocal, ALIGN_SIZE(rbmd::Real, 3),
                       H2D));

  CHECK_RUNTIME(MEMCPY(&_linked_cell_device_data_ptr->_d_cutoff, &this->_cutoff,
                       sizeof(rbmd::Real), H2D));
}

void LinkedCell::SortAtomsByCellKey() {
  // 使用 zip_iterator 组合多个数组
  auto zip_begin = thrust::make_zip_iterator(thrust::make_tuple(
      _device_data->_d_atoms_id.begin(),
      _device_data->_d_atoms_type.begin(),
      _device_data->_d_px.begin(),
      _device_data->_d_py.begin(),
      _device_data->_d_pz.begin(),
      _device_data->_d_vx.begin(),
      _device_data->_d_vy.begin(),
      _device_data->_d_vz.begin(),
      _device_data->_d_charge.begin()
  ));

  // 一次性排序所有数据
  thrust::stable_sort_by_key(_per_atom_cell_id.begin(),
                             _per_atom_cell_id.end(),
                             zip_begin);

  // 执行 MapAtomidToIdxOp
  op::MapAtomidToIdxOp<device::DEVICE_GPU> map_atomid_to_idx_op;
  map_atomid_to_idx_op(thrust::raw_pointer_cast(_atom_id_to_idx.data()),
                       raw_ptr(_device_data->_d_atoms_id), _total_atoms_num);
}

void LinkedCell::AllocDeviceMemory() {
  if (nullptr == this->_linked_cell_device_data_ptr) {
    CHECK_RUNTIME(MALLOC(&(this->_linked_cell_device_data_ptr),
                         sizeof(LinkedCellDeviceDataPtr)));
  }
}
