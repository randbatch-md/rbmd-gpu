#include "include/scheduler/memory_scheduler.h"

#include <thrust/copy.h>

#include "../data_manager/include/data_manager.h"

MemoryScheduler::MemoryScheduler()
    :  //_data_manager(DataManager::getInstance()),
      _device_data(DataManager::getInstance().getDeviceData()),
      _md_data(DataManager::getInstance().getMDData()),
      _structure_data(_md_data->_structure_data),
      _force_field_data(_md_data->_force_field_data),
      _structure_info_data(_md_data->_structure_info_data) {}

bool MemoryScheduler::asyncMemoryH2D() {
  auto& num_atoms = *(_structure_info_data->_num_atoms);
  auto& atoms_type = *(_structure_info_data->_num_atoms_type);
  auto& num_bonds = *(_structure_info_data->_num_bonds);
  auto& num_angles = *(_structure_info_data->_num_angles);
  auto& num_dihedrals = *(_structure_info_data->_num_dihedrals);
  auto& num_impropers = *(_structure_info_data->_num_impropers);

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

  auto& h_mass = _force_field_data->_h_mass;
  auto& h_atoms_id = _structure_data->_h_atoms_id;
  auto& h_atoms_type = _structure_data->_h_atoms_type;
  auto& h_molecular_id = _structure_data->_h_molecular_id;

  /// copy position
  _device_data->_d_px.resize(num_atoms);
  _device_data->_d_py.resize(num_atoms);
  _device_data->_d_pz.resize(num_atoms);
  _device_data->_d_shake_px.resize(num_atoms);
  _device_data->_d_shake_py.resize(num_atoms);
  _device_data->_d_shake_pz.resize(num_atoms);

  thrust::copy(h_px, h_px + num_atoms, _device_data->_d_px.begin());
  thrust::copy(h_py, h_py + num_atoms, _device_data->_d_py.begin());
  thrust::copy(h_pz, h_pz + num_atoms, _device_data->_d_pz.begin());

  /// copy flag
  _device_data->_d_flagX.resize(num_atoms,0.0);
  _device_data->_d_flagY.resize(num_atoms,0.0);
  _device_data->_d_flagZ.resize(num_atoms,0.0);
  _device_data->_d_unwarp_px.resize(num_atoms);
  _device_data->_d_unwarp_py.resize(num_atoms);
  _device_data->_d_unwarp_pz.resize(num_atoms);
  thrust::copy(h_px, h_px + num_atoms, _device_data->_d_unwarp_px.begin());
  thrust::copy(h_py, h_py + num_atoms, _device_data->_d_unwarp_py.begin());
  thrust::copy(h_pz, h_pz + num_atoms, _device_data->_d_unwarp_pz.begin());

  _device_data->_d_shake_vx.resize(num_atoms);
  _device_data->_d_shake_vy.resize(num_atoms);
  _device_data->_d_shake_vz.resize(num_atoms);
  // thrust::copy(h_flagX, h_flagX + num_atoms, _device_data->_d_flagX.begin());
  // thrust::copy(h_flagY, h_flagY + num_atoms, _device_data->_d_flagY.begin());
  // thrust::copy(h_flagZ, h_flagZ + num_atoms, _device_data->_d_flagZ.begin());

  /// copy velocity
  _device_data->_d_vx.resize(num_atoms);
  _device_data->_d_vy.resize(num_atoms);
  _device_data->_d_vz.resize(num_atoms);
  thrust::copy(h_vx, h_vx + num_atoms, _device_data->_d_vx.begin());
  thrust::copy(h_vy, h_vy + num_atoms, _device_data->_d_vy.begin());
  thrust::copy(h_vz, h_vz + num_atoms, _device_data->_d_vz.begin());



  // copy charge
  if (charge_data && charge_data->_h_charge)
  {
    _device_data->_d_charge.resize(num_atoms);
    thrust::copy(charge_data->_h_charge, charge_data->_h_charge + num_atoms,
      _device_data->_d_charge.begin());
  }
  else{_device_data->_d_charge.resize(num_atoms);}


  //resize force
  _device_data->_d_fx.resize(num_atoms);
  _device_data->_d_fy.resize(num_atoms);
  _device_data->_d_fz.resize(num_atoms);
  // thrust::copy(h_fx, h_fx + num_atoms, _device_data->_d_fx.begin());
  // thrust::copy(h_fy, h_fy + num_atoms, _device_data->_d_fy.begin());
  // thrust::copy(h_fz, h_fz + num_atoms, _device_data->_d_fz.begin());

  // _device_data->_d_prev_fx.resize(num_atoms);
  // _device_data->_d_prev_fy.resize(num_atoms);
  // _device_data->_d_prev_fz.resize(num_atoms);
  // _device_data->_d_pr_prev_fx.resize(num_atoms);
  // _device_data->_d_pr_prev_fy.resize(num_atoms);
  // _device_data->_d_pr_prev_fz.resize(num_atoms);


  _device_data->_d_force_ljcoul_x.resize(num_atoms);
  _device_data->_d_force_ljcoul_y.resize(num_atoms);
  _device_data->_d_force_ljcoul_z.resize(num_atoms);

  _device_data->_d_force_specialcoul_x.resize(num_atoms);
  _device_data->_d_force_specialcoul_y.resize(num_atoms);
  _device_data->_d_force_specialcoul_z.resize(num_atoms);


  _device_data->_d_force_kspace_x.resize(num_atoms);
  _device_data->_d_force_kspace_y.resize(num_atoms);
  _device_data->_d_force_kspace_z.resize(num_atoms);

  //bond
  _device_data->_d_force_bond_x.resize(num_atoms);
  _device_data->_d_force_bond_y.resize(num_atoms);
  _device_data->_d_force_bond_z.resize(num_atoms);

  _device_data-> _d_temp_atom_ids.resize(num_bonds * 2);
  _device_data->_d_temp_forces_bondx.resize(num_bonds * 2);
  _device_data->_d_temp_forces_bondy.resize(num_bonds * 2);
  _device_data->_d_temp_forces_bondz.resize(num_bonds * 2);
  //angle
  _device_data->_d_force_angle_x.resize(num_atoms);
  _device_data->_d_force_angle_y.resize(num_atoms);
  _device_data->_d_force_angle_z.resize(num_atoms);
  //dihedral
  _device_data->_d_force_dihedral_x.resize(num_atoms);
  _device_data->_d_force_dihedral_y.resize(num_atoms);
  _device_data->_d_force_dihedral_z.resize(num_atoms);
  //improper
  _device_data->_d_force_improper_x.resize(num_atoms);
  _device_data->_d_force_improper_y.resize(num_atoms);
  _device_data->_d_force_improper_z.resize(num_atoms);

  _device_data->_d_f2_x.resize(num_atoms);
  _device_data->_d_f2_y.resize(num_atoms);
  _device_data->_d_f2_z.resize(num_atoms);
  _device_data->_d_f3_x.resize(num_atoms);
  _device_data->_d_f3_y.resize(num_atoms);
  _device_data->_d_f3_z.resize(num_atoms);
  //copy virial
  _device_data->_d_flat_virial.resize(6*num_atoms);
  _device_data->_d_flat_virial_lj.resize(6*num_atoms);
  _device_data->_d_flat_virial_specialcoul.resize(6*num_atoms);
  _device_data->_d_flat_virial_kspace.resize(6*num_atoms);

  _device_data->_d_flat_virial_bond_list.resize(6*num_bonds);
  _device_data->_d_flat_virial_angle_list.resize(6*num_angles);
  _device_data->_d_flat_virial_dihedral_list.resize(6*num_dihedrals);
  _device_data->_d_flat_virial_improper_list.resize(6*num_impropers);

  _device_data->_d_flat_virial_bond_atom.resize(6*num_atoms);
  _device_data->_d_flat_virial_angle_atom.resize(6*num_atoms);
  _device_data->_d_flat_virial_dihedral_atom.resize(6*num_atoms);
  _device_data->_d_flat_virial_improper_atom.resize(6*num_atoms);

  _device_data->_d_virial.resize(6);
  _device_data->_d_virial_lj.resize(6);
  _device_data->_d_virial_specialcoul.resize(6);
  _device_data->_d_virial_kspace.resize(6);
  _device_data->_d_virial_bond.resize(6);
  _device_data->_d_virial_angle.resize(6);
  _device_data->_d_virial_dihedral.resize(6);
  _device_data->_d_virial_improper.resize(6);

  /// copy other
  _device_data->_d_mass.resize(atoms_type);
  _device_data->_d_atoms_id.resize(num_atoms);
  _device_data->_d_atoms_type.resize(num_atoms);
  //_device_data->_d_molecular_id.resize(num_atoms);

  thrust::copy(h_mass, h_mass + atoms_type,
             _device_data->_d_mass.begin());
  thrust::copy(h_atoms_id, h_atoms_id + num_atoms,
               _device_data->_d_atoms_id.begin());
  thrust::copy(h_atoms_type, h_atoms_type + num_atoms,
               _device_data->_d_atoms_type.begin());
  // thrust::copy(h_molecular_id, h_molecular_id + num_atoms,    //TODO
  //
  //              _device_data->_d_molecular_id.begin());

  // copy box

  //
  CHECK_RUNTIME(MALLOC(&_device_data->_d_erf_table, sizeof(ERFTable)));
  CHECK_RUNTIME(
      MEMCPY(_device_data->_d_erf_table, _md_data->_h_erf_table.get(), sizeof(ERFTable), H2D));
  return true;
}

bool MemoryScheduler::asyncMemoryD2H() { return true; }
