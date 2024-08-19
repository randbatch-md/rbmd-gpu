#include <thrust/copy.h>
#include "include/lj_memory_scheduler.h"
#include "include/data_manager.h"
#include "include/device_data.h"
#include "include/md_data.h"

bool LJMemoryScheduler::asyncMemoryH2D()
{
	auto& data_manager = DataManager::getInstance();
	auto& device_data = data_manager.getDeviceData();
	auto& md_data = data_manager.getMDData();
	auto& structure_data = md_data->_structure_data;
	auto& force_field = md_data->_force_field_data;
	auto& structure_info = md_data->_structure_info_data;

	auto& num_atoms = structure_info->_num_atoms;
	auto& num_atom_types = structure_info->_num_atoms_type;

	auto& h_px = structure_data->_h_px;
	auto& h_py = structure_data->_h_py;
	auto& h_pz = structure_data->_h_pz;

	auto& h_vx = structure_data->_h_vx;
	auto& h_vy = structure_data->_h_vy;
	auto& h_vz = structure_data->_h_vz;

	auto& h_fx = structure_data->_h_fx;
	auto& h_fy = structure_data->_h_fy;
	auto& h_fz = structure_data->_h_fz;

	auto& h_atoms_id = structure_data->_h_atoms_id;
	auto& h_atoms_type = structure_data->_h_atoms_type;
	auto& h_molecular_id = structure_data->_h_molecular_id;

	auto& h_eps = force_field->_h_eps;
	auto& h_mass = force_field->_h_mass;
	auto& h_sigma = force_field->_h_sigma;

	///cpoy position
	device_data->_d_px.resize(num_atoms);
	device_data->_d_py.resize(num_atoms);
	device_data->_d_pz.resize(num_atoms);
	thrust::copy(h_px, h_px + num_atoms, device_data->_d_px.begin());
	thrust::copy(h_py, h_py + num_atoms, device_data->_d_py.begin());
	thrust::copy(h_pz, h_pz + num_atoms, device_data->_d_pz.begin());

	///cpoy velocity
	device_data->_d_vx.resize(num_atoms);
	device_data->_d_vy.resize(num_atoms);
	device_data->_d_vz.resize(num_atoms);
	thrust::copy(h_vx, h_vx + num_atoms, device_data->_d_vx.begin());
	thrust::copy(h_vy, h_vy + num_atoms, device_data->_d_vy.begin());
	thrust::copy(h_vz, h_vz + num_atoms, device_data->_d_vz.begin());

	///cpoy force
	device_data->_d_fx.resize(num_atoms);
	device_data->_d_fy.resize(num_atoms);
	device_data->_d_fz.resize(num_atoms);
	thrust::copy(h_fx, h_fx + num_atoms, device_data->_d_fx.begin());
	thrust::copy(h_fy, h_fy + num_atoms, device_data->_d_fy.begin());
	thrust::copy(h_fz, h_fz + num_atoms, device_data->_d_fz.begin());

	///copy other
	device_data->_d_atoms_id.resize(num_atoms);
	device_data->_d_atoms_type.resize(num_atoms);
	device_data->_d_molecular_id.resize(num_atoms);
	thrust::copy(h_atoms_id, h_atoms_id + num_atoms, device_data->_d_atoms_id.begin());
	thrust::copy(h_atoms_type, h_atoms_type + num_atoms, device_data->_d_atoms_type.begin());
	thrust::copy(h_molecular_id, h_molecular_id + num_atoms, device_data->_d_molecular_id.begin());

	///copy force field
	device_data->_d_eps.resize(num_atom_types);
	device_data->_d_mass.resize(num_atom_types);
	device_data->_d_sigma.resize(num_atom_types);
	thrust::copy(h_eps, h_eps + num_atom_types, device_data->_d_eps.begin());
	thrust::copy(h_mass, h_mass + num_atom_types, device_data->_d_mass.begin());
	thrust::copy(h_sigma, h_sigma + num_atom_types, device_data->_d_sigma.begin());

	return true;
}

bool LJMemoryScheduler::asyncMemoryD2H()
{
	return true;
}
