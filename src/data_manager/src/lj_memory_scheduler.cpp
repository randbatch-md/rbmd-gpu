#include "include/scheduler/lj_memory_scheduler.h"

bool LJMemoryScheduler::asyncMemoryH2D() {
  if (false == MemoryScheduler::asyncMemoryH2D()) {
    // log
    return false;
  }

  auto& num_atoms_type = *(_structure_info_data->_num_atoms_type);
  auto& num_atoms = *(_structure_info_data-> _num_atoms);
  std::cout << "num_atoms_type" << num_atoms_type << std::endl;
  auto fd = std::dynamic_pointer_cast<LJForceFieldData>(_force_field_data);

  /// copy force field
  _device_data->_d_eps.resize(num_atoms_type);
  _device_data->_d_mass.resize(num_atoms_type);
  _device_data->_d_sigma.resize(num_atoms_type);

  //_device_data->_d_charge.resize(num_atoms);

  /// eps
  thrust::copy(fd->_h_eps, fd->_h_eps + num_atoms_type,
               _device_data->_d_eps.begin());

  /// mass
  thrust::copy(fd->_h_mass, fd->_h_mass + num_atoms_type,
               _device_data->_d_mass.begin());

  /// sigma
  thrust::copy(fd->_h_sigma, fd->_h_sigma + num_atoms_type,
               _device_data->_d_sigma.begin());

  /// charge
  // thrust::copy(fd->_h_charge, fd->_h_charge + num_atoms,
  // _device_data->_d_charge.begin());

  return true;
}

bool LJMemoryScheduler::asyncMemoryD2H() { return true; }
