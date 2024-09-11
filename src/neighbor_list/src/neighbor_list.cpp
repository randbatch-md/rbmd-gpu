#include "neighbor_list/include/neighbor_list.h"

#include "common/rbmd_define.h"

NeighborList::NeighborList(rbmd::Id total_atoms_num, bool is_half) {
  CHECK_RUNTIME(
      MALLOCHOST(reinterpret_cast<void **>(&(this->_h_is_half)), sizeof(bool)));
  *(this->_h_is_half) = is_half;
  this->_d_neighbor_num.resize(total_atoms_num);
  this->_d_max_neighbor_num.resize(total_atoms_num);
  this->_start_idx.resize(total_atoms_num);
  this->_end_idx.resize(total_atoms_num);
  CHECK_RUNTIME(MALLOC(&this->_d_is_half, sizeof(bool)));
  CHECK_RUNTIME(MEMCPY(this->_d_is_half, this->_h_is_half, sizeof(bool), H2D));
}
