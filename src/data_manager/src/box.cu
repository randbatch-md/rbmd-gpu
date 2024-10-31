#include "model/box.h"

void Box::Init(BoxType box_type, const rbmd::Real coord_min[3],
               const rbmd::Real coord_max[3], const bool pbc[3]) {
  this->_type = box_type;
  for (int i = 0; i < 3; ++i) {
    this->_coord_min[i] = coord_min[i];
    this->_coord_max[i] = coord_max[i];
    this->_length[i] = coord_max[i] - coord_min[i];
    this->_length_inv[i] = 1 / _length[i];
  }
  this->_length[3] = 0;  // yz
  this->_length[4] = 0;  // xz
  this->_length[5] = 0;  // xy

  this->_pbc_x = pbc[0];
  this->_pbc_y = pbc[1];
  this->_pbc_z = pbc[2];
}