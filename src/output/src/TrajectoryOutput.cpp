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

TrajectoryOutput::TrajectoryOutput()
: _interval(DataManager::getInstance().getConfigData()->Get<rbmd::Id>("interval", "outputs", "trajectory_out"))
{
    spdlog::init_thread_pool(10000, 1);
    auto async_file_logger = spdlog::basic_logger_mt<spdlog::async_factory>("async_file_logger", "rbmd.trj");
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
            std::vector<rbmd::Real> h_px(_num_atoms);
            std::vector<rbmd::Real> h_py(_num_atoms);
            std::vector<rbmd::Real> h_pz(_num_atoms);
            std::vector<rbmd::Real> h_atoms_type(_num_atoms);
            rbmd::Range h_rang;

            thrust::copy(_device_data->_d_px.begin(), _device_data->_d_px.end(), h_px.begin());
            thrust::copy(_device_data->_d_py.begin(), _device_data->_d_py.end(), h_py.begin());
            thrust::copy(_device_data->_d_pz.begin(), _device_data->_d_pz.end(), h_pz.begin());
            thrust::copy(_device_data->_d_atoms_type.begin(), _device_data->_d_atoms_type.end(), h_atoms_type.begin());

            CHECK_RUNTIME(MEMCPY(&h_rang, _structure_info_data->_range, sizeof(rbmd::Range), D2H));

            std::cout<<"h_rang value: "<< h_rang[0][0] << " " << h_rang[0][1] << std::endl
                     << h_rang[1][0] << " " << h_rang[1][1] << std::endl
                     << h_rang[2][0] << " " << h_rang[2][1] << std::endl;

            spdlog::info("ITEM: TIMESTEP");
            spdlog::info("{}", test_current_step);
            spdlog::info("ITEM: NUMBER OF ATOMS");
            spdlog::info("{}", _num_atoms);
            spdlog::info("ITEM: BOX BOUNDS pp pp pp");
            spdlog::info("{} {}", h_rang[0][0], h_rang[0][1]);
            spdlog::info("{} {}", h_rang[1][0], h_rang[1][1]);
            spdlog::info("{} {}", h_rang[2][0], h_rang[2][1]);
            spdlog::info("ITEM: ATOMS id type x y z");

            for (auto i = 0; i < _num_atoms; ++i)
            {
                spdlog::info("{} {} {} {} {}", i + 1, h_atoms_type[i] + 1, h_px[i], h_py[i], h_pz[i]);
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
    if (test_current_step < 1)
    {
        return false;
    }
    return test_current_step % _interval == 0;
}