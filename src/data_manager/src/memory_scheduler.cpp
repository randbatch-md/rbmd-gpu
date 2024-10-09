#include "include/scheduler/memory_scheduler.h"
#include "../data_manager/include/data_manager.h"

#include <thrust/copy.h>

MemoryScheduler::MemoryScheduler()
    : //_data_manager(DataManager::getInstance()),
      _device_data(DataManager::getInstance().getDeviceData()),
      _md_data(DataManager::getInstance().getMDData()),
      _structure_data(_md_data->_structure_data),
      _force_field_data(_md_data->_force_field_data),
      _structure_info_data(_md_data->_structure_info_data) {}

bool MemoryScheduler::asyncMemoryH2D() {
  auto& num_atoms = *(_structure_info_data->_num_atoms);

  auto& h_px = _structure_data->_h_px;
  auto& h_py = _structure_data->_h_py;
  auto& h_pz = _structure_data->_h_pz;
  //
  auto& h_flagX = _structure_data->_h_flagX;
  auto& h_flagY = _structure_data->_h_flagY;
  auto& h_flagZ = _structure_data->_h_flagZ;

  auto& h_vx = _structure_data->_h_vx;
  auto& h_vy = _structure_data->_h_vy;
  auto& h_vz = _structure_data->_h_vz;

  auto& h_fx = _structure_data->_h_fx;
  auto& h_fy = _structure_data->_h_fy;
  auto& h_fz = _structure_data->_h_fz;

  //
  ChargeStructureData* charge_data = dynamic_cast<ChargeStructureData*>(_structure_data.get());
  auto& h_charge = charge_data->_h_charge;

  auto& h_evdwl= _structure_data->_h_evdwl;

  auto& h_atoms_id = _structure_data->_h_atoms_id;
  auto& h_atoms_type = _structure_data->_h_atoms_type;
  auto& h_molecular_id = _structure_data->_h_molecular_id;

  /// cpoy position
  _device_data->_d_px.resize(num_atoms);
  _device_data->_d_py.resize(num_atoms);
  _device_data->_d_pz.resize(num_atoms);
  thrust::copy(h_px, h_px + num_atoms, _device_data->_d_px.begin());
  thrust::copy(h_py, h_py + num_atoms, _device_data->_d_py.begin());
  thrust::copy(h_pz, h_pz + num_atoms, _device_data->_d_pz.begin());

  /// cpoy flag
  _device_data->_d_flagX.resize(num_atoms);
  _device_data->_d_flagY.resize(num_atoms);
  _device_data->_d_flagZ.resize(num_atoms);
  // thrust::copy(h_flagX, h_flagX + num_atoms, _device_data->_d_flagX.begin());
  // thrust::copy(h_flagY, h_flagY + num_atoms, _device_data->_d_flagY.begin());
  // thrust::copy(h_flagZ, h_flagZ + num_atoms, _device_data->_d_flagZ.begin());


  /// cpoy velocity
  _device_data->_d_vx.resize(num_atoms);
  _device_data->_d_vy.resize(num_atoms);
  _device_data->_d_vz.resize(num_atoms);
  thrust::copy(h_vx, h_vx + num_atoms, _device_data->_d_vx.begin());
  thrust::copy(h_vy, h_vy + num_atoms, _device_data->_d_vy.begin());
  thrust::copy(h_vz, h_vz + num_atoms, _device_data->_d_vz.begin());



  // cpoy charge
  _device_data->_d_charge.resize(num_atoms);
  thrust::copy(h_charge, h_charge + num_atoms, _device_data->_d_charge.begin());

  /// cpoy force
  _device_data->_d_fx.resize(num_atoms);
  _device_data->_d_fy.resize(num_atoms);
  _device_data->_d_fz.resize(num_atoms);
   // thrust::copy(h_fx, h_fx + num_atoms, _device_data->_d_fx.begin());
   // thrust::copy(h_fy, h_fy + num_atoms, _device_data->_d_fy.begin());
   // thrust::copy(h_fz, h_fz + num_atoms, _device_data->_d_fz.begin());

  _device_data->_d_force_ljcoul_x.resize(num_atoms);
  _device_data->_d_force_ljcoul_y.resize(num_atoms);
  _device_data->_d_force_ljcoul_z.resize(num_atoms);

  _device_data->_d_force_ewald_x.resize(num_atoms);
  _device_data->_d_force_ewald_y.resize(num_atoms);
  _device_data->_d_force_ewald_z.resize(num_atoms);

  _device_data->_d_virial_xx.resize(num_atoms);
  _device_data->_d_virial_yy.resize(num_atoms);
  _device_data->_d_virial_zz.resize(num_atoms);
  _device_data->_d_virial_xy.resize(num_atoms);
  _device_data->_d_virial_xy.resize(num_atoms);
  _device_data->_d_virial_yz.resize(num_atoms);


   //cpoy evdwl
  _device_data->_d_evdwl.resize(num_atoms);
  //thrust::copy(h_evdwl, h_evdwl + num_atoms, _device_data->_d_evdwl.begin());

  /// copy other
  _device_data->_d_atoms_id.resize(num_atoms);
  _device_data->_d_atoms_type.resize(num_atoms);
  _device_data->_d_molecular_id.resize(num_atoms);
  thrust::copy(h_atoms_id, h_atoms_id + num_atoms,
               _device_data->_d_atoms_id.begin());
  thrust::copy(h_atoms_type, h_atoms_type + num_atoms,
               _device_data->_d_atoms_type.begin());
  // thrust::copy(h_molecular_id, h_molecular_id + num_atoms,    //TODO 有好多h没有，只需resize
  //              _device_data->_d_molecular_id.begin());

  // copy box
  CHECK_RUNTIME(MALLOC(&_device_data->_d_box, sizeof(Box)));
  CHECK_RUNTIME(MEMCPY(_device_data->_d_box,_md_data->_h_box.get() , sizeof(Box), H2D));
  return true;
}

bool MemoryScheduler::asyncMemoryD2H() { return true; }
