#include "ljforce.h"
#include "ljforce_op/ljforce_op.h"
#include <thrust/device_vector.h>

LJForce::LJForce()
    : list(NeighborList(_structure_info_data->_num_atoms)) {};

void LJForce::Init()
{
    _num_atoms = _structure_info_data->_num_atoms;
}

void LJForce::Execute()
{
    LJForce::Init();
    //

    //BuildNeighbor(num_atoms, box, start_id, end_id, id_verletlist);

	op::LJforceOp<device::DEVICE_GPU> LJforceOp;

    LJforceOp(box, _num_atoms,
        thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
        thrust::raw_pointer_cast(_device_data->_d_molecular_id.data()),
        thrust::raw_pointer_cast(_device_data->_d_eps.data()),
        thrust::raw_pointer_cast(_device_data->_d_sigma.data()),
        thrust::raw_pointer_cast(list._start_idx.data()),
        thrust::raw_pointer_cast(list._end_idx.data()),
        thrust::raw_pointer_cast(list._d_neighbor_num.data()),
        thrust::raw_pointer_cast(_device_data->_d_px.data()),
        thrust::raw_pointer_cast(_device_data->_d_py.data()),
        thrust::raw_pointer_cast(_device_data->_d_pz.data()),
        thrust::raw_pointer_cast(_device_data->_d_fx.data()),
        thrust::raw_pointer_cast(_device_data->_d_fy.data()),
        thrust::raw_pointer_cast(_device_data->_d_fz.data()),
        thrust::raw_pointer_cast(_device_data->_d_evdwl.data()));

}

