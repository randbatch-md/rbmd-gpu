#include "TrajectoryOutput.h"
#include <thrust/device_ptr.h>
#include "data_manager.h"
#include "device_types.h"
#include "neighbor_list/include/linked_cell/linked_cell_locator.h"
#include "unit_factor.h"
#include "output_op.h"
#include <thrust/copy.h>
#include "spdlog/spdlog.h"
#include <spdlog/async.h>
#include <spdlog/sinks/basic_file_sink.h>
#include "box.h"
#include "simulate.h"
//#include <filesystem>
#include <fstream>
TrajectoryOutput::TrajectoryOutput()
: _interval(DataManager::getInstance().getConfigData()->Get<rbmd::Id>
  ("interval", "outputs", "trajectory_out"))
{
    //std::filesystem::remove("rbmd.trj");
    std::remove("rbmd.trj");
    spdlog::init_thread_pool(10000, 1);
    auto async_file_logger = spdlog::basic_logger_mt
      <spdlog::async_factory>("async_file_logger", "rbmd.trj", false);
    async_file_logger->set_level(spdlog::level::debug);
    async_file_logger->set_pattern("%v");
    spdlog::set_default_logger(async_file_logger);
};

void TrajectoryOutput::Init() 
{
    _num_atoms = *(_structure_info_data->_num_atoms);
}

void TrajectoryOutput::Execute()
{
    if (ShouldOutput())
    {
        try
        {
            auto atom_id_to_idx =
              LinkedCellLocator::GetInstance().GetLinkedCell()->_atom_id_to_idx;
            //xyz
            std::vector<rbmd::Real> h_px(_num_atoms);
            std::vector<rbmd::Real> h_py(_num_atoms);
            std::vector<rbmd::Real> h_pz(_num_atoms);
            std::vector<rbmd::Real> h_atoms_type(_num_atoms);
            //vx vy vz
            thrust::host_vector<rbmd::Real> h_vx(_device_data->_d_vx);
            thrust::host_vector<rbmd::Real> h_vy(_device_data->_d_vy);
            thrust::host_vector<rbmd::Real> h_vz(_device_data->_d_vz);

            //charge
            thrust::host_vector<rbmd::Real> h_charge(_device_data->_d_charge);

            auto box =  DataManager::getInstance().getMDData()->_box;

            thrust::copy(_device_data->_d_px.begin(), _device_data->_d_px.end(), h_px.begin());
            thrust::copy(_device_data->_d_py.begin(), _device_data->_d_py.end(), h_py.begin());
            thrust::copy(_device_data->_d_pz.begin(), _device_data->_d_pz.end(), h_pz.begin());
            thrust::copy(_device_data->_d_atoms_type.begin(), _device_data->_d_atoms_type.end(), h_atoms_type.begin());

            spdlog::info("ITEM: TIMESTEP");
            spdlog::info("{}", test_current_step);
            spdlog::info("ITEM: NUMBER OF ATOMS");
            spdlog::info("{}", _num_atoms);
            spdlog::info("ITEM: BOX BOUNDS pp pp pp");
            spdlog::info("{} {}", box->_coord_min[0], box->_coord_max[0]);
            spdlog::info("{} {}", box->_coord_min[1], box->_coord_max[1]);
            spdlog::info("{} {}", box->_coord_min[2], box->_coord_max[2]);
            spdlog::info("ITEM: ATOMS id type q x y z vx vy vz");

            for (auto i = 0; i < _num_atoms; ++i)
            {
                auto idx = atom_id_to_idx[i];
                spdlog::info("{} {} {} {} {} {} {} {} {}" ,
                  i + 1, h_atoms_type[idx] + 1, h_charge[idx],
                  h_px[idx], h_py[idx],h_pz[idx],
                  h_vx[idx],h_vy[idx],h_vz[idx]);
            }
        }
        catch (const std::exception& e)
        {
            throw std::runtime_error("TrajectoryOutput is error ");
        }
    }
}

bool TrajectoryOutput::ShouldOutput()
{
    // include the original date of step = 0;
    //if (test_current_step < 1)
    //{
    //    return false;
    //}
    return test_current_step % _interval == 0;
}