#include "mace_memory_scheduler.h"



#include "data_manager.h"
#include "mace_force_field_data.h"

bool MACEMemoryScheduler::asyncMemoryH2D() {
  if (false == MemoryScheduler::asyncMemoryH2D()) {
    // log
    return false;
  }

  auto& num_atoms_type = *(_structure_info_data->_num_atoms_type);
  auto& num_atoms = *(_structure_info_data->_num_atoms);
  //auto sd = std::dynamic_pointer_cast<FullStructureData>(_structure_data);
  auto fd = std::dynamic_pointer_cast<ForceFieldData>(_force_field_data);


  /// mass
  //   _device_data->_d_mass.resize(num_atoms_type);
  // std::cout << " size"<<_device_data->_d_mass.size() << std::endl;
  // thrust::copy(fd->_h_mass, fd->_h_mass + num_atoms_type,
  //              _device_data->_d_mass.begin());
  // std::cout << " kkkkk"<<num_atoms << std::endl;
  return true;
}

bool MACEMemoryScheduler::asyncMemoryD2H() { return true; }
