#include "./include/scheduler/eam_memory_scheduler.h"

bool EAMMemoryScheduler::asyncMemoryH2D() {
  if (false == MemoryScheduler::asyncMemoryH2D()) {
    // log
    return false;
  }
  auto fd = std::dynamic_pointer_cast<EAMForceFieldData>(_force_field_data);
  auto& num_atoms = *(_structure_info_data->_num_atoms);

  return true;
}

bool EAMMemoryScheduler::asyncMemoryD2H() { return true; }
