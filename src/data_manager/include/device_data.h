#pragma once
#include <thrust/device_vector.h>
#include "types.h"

class DeviceData
{
public:
	rbmd::Real test;
	thrust::device_vector<rbmd::Real> _d_px;
	thrust::device_vector<rbmd::Real> _d_py;
	thrust::device_vector<rbmd::Real> _d_pz;

	thrust::device_vector<rbmd::Id> _d_atoms_id;
	thrust::device_vector<rbmd::Id> _d_atoms_type;
	thrust::device_vector<rbmd::Id> _d_molecular_id;

	thrust::device_vector<rbmd::Real> _d_vx;
	thrust::device_vector<rbmd::Real> _d_vy;
	thrust::device_vector<rbmd::Real> _d_vz;

	thrust::device_vector<rbmd::Real> _d_fx;
	thrust::device_vector<rbmd::Real> _d_fy;
	thrust::device_vector<rbmd::Real> _d_fz;

	thrust::device_vector<rbmd::Real> _d_mass;
	thrust::device_vector<rbmd::Real> _d_eps;
	thrust::device_vector<rbmd::Real> _d_sigma;
};
