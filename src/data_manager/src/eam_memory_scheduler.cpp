#include "./include/scheduler/eam_memory_scheduler.h"

bool EAMMemoryScheduler::asyncMemoryH2D() {
  if (false == MemoryScheduler::asyncMemoryH2D()) {
    // log
    return false;
  }
  auto fd = std::dynamic_pointer_cast<EAMForceFieldData>(_force_field_data);

  auto& num_atoms = *(_structure_info_data->_num_atoms);
  /// copy force field
  // _device_data->_d_frho.resize(fd->_nrho);
  // _device_data->_d_rhor.resize(fd->_nrho);
  // _device_data->_d_z2r.resize(fd->_nr);

  // // /// frho
  // thrust::copy(fd->_h_frho.begin(), fd->_h_frho.end(),
  //              _device_data->_d_frho.begin());
  //
  // /// rhor
  // thrust::copy(fd->_h_rhor.begin(), fd->_h_rhor.end(),
  //              _device_data->_d_rhor.begin());
  //
  // /// z2r
  // thrust::copy(fd->_h_zr.begin(), fd->_h_zr.end(),
  //              _device_data->_d_zr.begin());

  _device_data->_eam_rho.resize(num_atoms);
  _device_data->_eam_fp.resize(num_atoms);
  return true;
}

bool EAMMemoryScheduler::asyncMemoryD2H() { return true; }
