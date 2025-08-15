#include "nvt_ensemble.h"

#include <chrono>  //

#include "data_manager.h"
#include "model/md_data.h"

#include "default_position_controller.h"
#include "default_velocity_controller.h"
#include "berendsen_controller.h"
#include "langevin_controller.h"
#include "rescale_controller.h"
#include "nose_hoover_controller.h"
#include "cvff.h"
#include "lj_cut_coul_kspace.h"
#include "lj.h"
#include "tersoff.h"
#include "eam.h"
#include "neighbor_list/include/linked_cell/linked_cell_locator.h"
#include "shake_controller.h"
#include "output/include/Logger.hpp"
NVTensemble::NVTensemble()
{
  _position_controller = std::make_shared<DefaultPositionController>();
  _velocity_controller = std::make_shared<DefaultVelocityController>();

  // Unified  Force Field Controller
  static const std::unordered_map<std::string, std::function<std::shared_ptr<Force>()>>
  force_map = {
    {"CVFF", [&]() { return std::make_shared<CVFF>(); }},
    {"LJ/CUT", [&]() { return std::make_shared<LJ>(); }},
    {"LJ/CUT/COUL/LONG", [&]() { return std::make_shared<LJCutCoulKspace>(); }},
    {"EAM", [&]() { return std::make_shared<EAM>(); }},
    {"Tersoff", [&]() { return std::make_shared<TerSoff>(); }}
  };

  //force_type
  auto force_type = DataManager::getInstance().getConfigData()->Get<std::string>
  ("type", "hyper_parameters", "force_field");
  if (auto it = force_map.find(force_type); it != force_map.end())
  {
    _force_controller = it->second();
  }
  else {
    Logger::Instance().error("Unsupported force field type: {}", force_type);
  }

  // unified temperature controller
  _temp_ctrl_type = DataManager::getInstance().getConfigData()->Get
  <std::string>("temp_ctrl_type", "execution");

  if ("RESCALE" == _temp_ctrl_type) {
    _temperature_controller = std::make_shared<RescaleController>();
  }
  else if ("BERENDSEN" == _temp_ctrl_type) {
    _temperature_controller = std::make_shared<BerendsenController>();
  }
  else if ("LANGEVIN" == _temp_ctrl_type) {
    _temperature_controller = std::make_shared<LangevinController>();
  }
  else if ("NOSE_HOOVER" == _temp_ctrl_type) {
    _NoseHoover_controller = std::make_shared<NoseHooverController>();
  }
  else {
    Logger::Instance().error("\033[31m Unsupported temp_ctrl_type: {}\033[0m", _temp_ctrl_type );
    exit(EXIT_FAILURE); //
  }

  // //
  // _NoseHoover_controller = std::make_shared<NoseHooverController>();

  //shake
  _shake_controller = std::make_shared<ShakeController>();
  _integration_type = DataManager::getInstance().getConfigData()->Get
<std::string>("integration_type", "execution");

}

void NVTensemble::Init() {
  _position_controller->Init();
  _velocity_controller->Init();

  _force_controller->Init();
  _force_controller->Execute();
  _shake_controller->Init();

  if (_temperature_controller) {
    _temperature_controller->Init();
  }

  if(_NoseHoover_controller) {
    _NoseHoover_controller->Init();
  }
}

void NVTensemble::Presolve() {}

void NVTensemble::Solve() {
   bool use_shake = DataManager::getInstance().getConfigData()->GetJudge
    <bool>("fix_shake", "hyper_parameters", "extend");
   if ("vv" ==_integration_type ) {
     if("NOSE_HOOVER" == _temp_ctrl_type)
     {
       _NoseHoover_controller->InitialIntegrate();//_velocity_controller->Update();
       //_position_controller->Update();
       if (true == use_shake)
       {
         _shake_controller->ShakeA();
       }

       _force_controller->Execute();

       _NoseHoover_controller->FinalIntegrate(); //_velocity_controller->Update();

       if (true == use_shake)
       {
         _shake_controller->ShakeB();
       }
     }
     else
     {
       auto start = std::chrono::high_resolution_clock::now();

   _velocity_controller->Update();

   _position_controller->Update();

   bool use_shake = DataManager::getInstance().getConfigData()->GetJudge<bool>
   ( "fix_shake", "hyper_parameters", "extend");; //TODO: json file
   if (use_shake)
   {
     _shake_controller->ShakeA();
   }

       _force_controller->Execute();

       if ("LANGEVIN"==DataManager::getInstance().getConfigData()->Get<std::string>
         ("temp_ctrl_type", "execution"))
       {
         _temperature_controller->Update();
       }

       _velocity_controller->Update();

       if (use_shake)
       {
         _shake_controller->ShakeB();
       }

       _temperature_controller->ComputeTemperature();

       if ("LANGEVIN" == DataManager::getInstance().getConfigData()->Get<std::string>
         ("temp_ctrl_type", "execution"))
         return;

       _temperature_controller->Update();

       CHECK_RUNTIME(DEVICESYNC());
       auto end = std::chrono::high_resolution_clock::now();
       std::chrono::duration<rbmd::Real> duration = end - start;
       _total_time_vv += duration.count();
       _total_steps++;
       std::cout << "VV单步耗时: " << duration.count() << "秒 | 平均耗时: " << (
         _total_time_vv / _total_steps) << "秒" << std::endl;
     }
   }
   else if ("test" ==_integration_type) {
     auto start = std::chrono::high_resolution_clock::now();
     _velocity_controller->Update1();
     _position_controller->Update1();
     _force_controller->Execute();
     _velocity_controller->Update2();
     _position_controller->Update2();
     _force_controller->Execute();
     _velocity_controller->Update3();
     _position_controller->Update3();
     _force_controller->Execute();
     _velocity_controller->Update4();
     _position_controller->Update4();
     _force_controller->Execute();

     _temperature_controller->ComputeTemperature();
     _temperature_controller->Update();
     CHECK_RUNTIME(DEVICESYNC());
     auto end = std::chrono::high_resolution_clock::now();
     std::chrono::duration<rbmd::Real> duration = end - start;
     _total_time_prk3c += duration.count();
     _total_steps++;
     std::cout << "PRK3C单步耗时: " << duration.count() << "秒 | 平均耗时: " << (
       _total_time_prk3c / _total_steps) << "秒" << std::endl;
    }
   else if ("vl" ==_integration_type) {
     // auto _device_data = DataManager::getInstance().getDeviceData();
     _position_controller->Updatevl();
     // thrust::host_vector<rbmd::Real> h_vx = _device_data->_d_vx;
     // printf("Before Execute: h_vx[1] = %f\n", h_vx[1]);
     // thrust::host_vector<rbmd::Real> h_fx = _device_data->_d_fx;
     // printf("Before Execute: h_fx[1] = %f\n", h_fx[1]);
     // thrust::host_vector<rbmd::Real> h_px = _device_data->_d_px;
     // printf("Before Execute: h_px[1] = %f\n", h_px[1]);
     // thrust::host_vector<rbmd::Real> h_px_prev = _device_data->_d_prev_px;
     // printf("Before Execute: h_px_prev[1] = %f\n", h_px_prev[1]);
     _force_controller->Execute();
     // thrust::host_vector<rbmd::Real> g_fx = _device_data->_d_fx;
     // printf("After Execute: g_fx[1] = %f\n", g_fx[1]);
     // thrust::host_vector<rbmd::Real> g_px = _device_data->_d_px;
     // printf("After Execute: g_px[1] = %f\n", g_px[1]);
     // thrust::host_vector<rbmd::Real> g_px_prev = _device_data->_d_prev_px;
     // printf("After Execute: g_px_prev[1] = %f\n", g_px_prev[1]);
     _velocity_controller->Updatevl();
     // thrust::host_vector<rbmd::Real> g_vx = _device_data->_d_vx;
     // printf("After Execute: g_vx[1] = %f\n", g_vx[1]);


     _temperature_controller->ComputeTemperature();
     _temperature_controller->Update();
   }
   else if ("bm" ==_integration_type) {
     // auto _device_data = DataManager::getInstance().getDeviceData();
     // auto start = std::chrono::high_resolution_clock::now();
     _position_controller->Updatebm();
   //   thrust::host_vector<rbmd::Real> h_vx = _device_data->_d_vx;
   //   printf("Before Execute: h_vx[1] = %f\n", h_vx[1]);
   //   thrust::host_vector<rbmd::Real> h_fx = _device_data->_d_fx;
   //   printf("Before Execute: h_fx[1] = %f\n", h_fx[1]);
   //   thrust::host_vector<rbmd::Real> h_px = _device_data->_d_px;
   //   printf("Before Execute: h_px[1] = %f\n", h_px[1]);
   //   thrust::host_vector<rbmd::Real> h_fx_prev = _device_data->_d_prev_fx;
   //   printf("Before Execute: h_fx_prev[1] = %f\n", h_fx_prev[1]);
   //
   //   auto atom_id_to_idx =
   //     LinkedCellLocator::GetInstance().GetLinkedCell()->_atom_id_to_idx;
   //   auto num_atoms = 8000;
   //   thrust::host_vector<rbmd::Real> h_f_x(num_atoms);
   //   thrust::host_vector<rbmd::Real> h_prev_fx(num_atoms);
   //
   //   thrust::copy(_device_data->_d_fx.begin(),
   // _device_data->_d_fx.end(), h_f_x.begin());
   //   thrust::copy(_device_data->_d_prev_fx.begin(),
   //     _device_data->_d_prev_fx.end(), h_prev_fx.begin());
   //
   //   std::ofstream output_file1("output_force1.txt");
   //   for (size_t i = 0; i < h_f_x.size(); ++i)
   //   {
   //     auto id = atom_id_to_idx[i];
   //     output_file1 << i << " " << h_f_x[id] << " "<<h_prev_fx[id] << std::endl;
   //   }
   //   output_file1.close();


     _force_controller->Execute();

   //   thrust::host_vector<rbmd::Real> h_fx2(num_atoms);
   //   thrust::host_vector<rbmd::Real> h_prev_fx2(num_atoms);
   //
   //
   //   thrust::copy(_device_data->_d_fx.begin(),
   // _device_data->_d_fx.end(), h_fx2.begin());
   //   thrust::copy(_device_data->_d_prev_fx.begin(),
   //     _device_data->_d_prev_fx.end(), h_prev_fx2.begin());
   //
   //
   //   std::ofstream output_file2("output_force2.txt");
   //   for (size_t i = 0; i < h_fx2.size(); ++i)
   //   {
   //     auto id = atom_id_to_idx[i];
   //     output_file2 << i << " " << h_fx2[id] << " "<<h_prev_fx2[id] << std::endl;
   //
   //   }
   //   output_file2.close();
   //
   //
   //
   //   thrust::host_vector<rbmd::Real> g_fx = _device_data->_d_fx;
   //   printf("After Execute: g_fx[1] = %f\n", g_fx[1]);
   //   thrust::host_vector<rbmd::Real> g_px = _device_data->_d_px;
   //   printf("After Execute: g_px[1] = %f\n", g_px[1]);
   //   thrust::host_vector<rbmd::Real> g_fx_prev = _device_data->_d_prev_fx;
   //   printf("After Execute: h_fx_prev[1] = %f\n", g_fx_prev[1]);
     _velocity_controller->Updatebm();
     // thrust::host_vector<rbmd::Real> g_vx = _device_data->_d_vx;
     // printf("After Execute: g_vx[1] = %f\n", g_vx[1]);
     _temperature_controller->ComputeTemperature();
     _temperature_controller->Update();
     // CHECK_RUNTIME(DEVICESYNC());
     // auto end = std::chrono::high_resolution_clock::now();
     // std::chrono::duration<rbmd::Real> duration = end - start;
     // _total_time_beeman += duration.count();
     // _total_steps++;
     // std::cout << "Beeman单步耗时: " << duration.count() << "秒 | 平均耗时: " << (
     //   _total_time_beeman / _total_steps) << "秒" << std::endl;
   }
}

void NVTensemble::Postsolve() {
  // if (_total_steps > 0) {
  //   Logger::Instance().info("===== 平均耗时统计 =====");
  //   if (_total_time_vv > 0) {
  //     Logger::Instance().info("VV算法平均耗时: {} 秒/步", _total_time_vv / _total_steps);
  //   }
  //   if (_total_time_beeman > 0) {
  //     Logger::Instance().info("Beeman算法平均耗时: {} 秒/步", _total_time_beeman / _total_steps);
  //   }
  //   if (_total_time_prk3c > 0) {
  //     Logger::Instance().info("PRK3C算法平均耗时: {} 秒/步", _total_time_prk3c / _total_steps);
  //   }
  // }
}
