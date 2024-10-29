#include "nvt_ensemble.h"

#include <chrono>  // 添加计时功能的库

#include "default_position_controller.h"
#include "default_velocity_controller.h"
#include "ljforce.h"
#include "lj_cut_coul_kspace_force.h"
#include "cvff.h"
#include "rescale_controller.h"
#include "berendsen_controller.h"
#include "nose_hoover_controller.h"
#include "shake_controller.h"
#define DEBUG
#include "data_manager.h"
#include "model/md_data.h"
#include "data_manager.h"
NVTensemble::NVTensemble()
{
  _position_controller = std::make_shared<DefaultPositionController>();
  _velocity_controller = std::make_shared<DefaultVelocityController>();
  _force_controller = std::make_shared<CVFF>(); // TODO: json file forcetype
  _temperature_controller = std::make_shared<BerendsenController>();
  _shake_controller = std::make_shared<ShakeController>();
}

void NVTensemble::Init() {
  _position_controller->Init();
  _velocity_controller->Init();
  _temperature_controller->Init();

  _force_controller->Init();
  _force_controller->Execute();
}

void NVTensemble::Presolve() {}

void NVTensemble::Solve() {
  auto start = std::chrono::high_resolution_clock::now();

  _velocity_controller->Update();

#ifdef DEBUG
  std::cout<<"---------第一次更新速度--------------"<<std::endl;
  std::shared_ptr<StructureInfoData> _structure_info_data;
  auto _device_data=DataManager::getInstance().getDeviceData();
  _structure_info_data= DataManager::getInstance().getMDData()->_structure_info_data;
  std::vector<rbmd::Real> h_vx(300);
  std::vector<rbmd::Real> h_px(300);
  thrust::copy(_device_data->_d_vx.begin(), _device_data->_d_vx.end(), h_vx.begin());
  thrust::copy(_device_data->_d_px.begin(), _device_data->_d_px.end(), h_px.begin());

  for (int j = 0; j < 10; ++j) {std::cout<<h_vx[j]<<" , ";}
  std::cout<<std::endl;
  for (int j = 0; j < 10; ++j) {std::cout<<h_px[j]<<" , ";}
  std::cout<<std::endl;
#endif

  _position_controller->Update();

#ifdef DEBUG
    std::cout<<"---------更新位置--------------"<<std::endl;
    _structure_info_data= DataManager::getInstance().getMDData()->_structure_info_data;
    thrust::copy(_device_data->_d_vx.begin(), _device_data->_d_vx.end(), h_vx.begin());
    thrust::copy(_device_data->_d_px.begin(), _device_data->_d_px.end(), h_px.begin());

    for (int j = 0; j < 10; ++j) {std::cout<<h_vx[j]<<" , ";}
    std::cout<<std::endl;
    for (int j = 0; j < 10; ++j) {std::cout<<h_px[j]<<" , ";}
    std::cout<<std::endl;
#endif

  bool use_shake = DataManager::getInstance().getConfigData()->GetJudge<bool>( "fix_shake", "hyper_parameters", "extend");; //TODO: json file
  if (use_shake)
  {
    _shake_controller->ShakeA();
  }

#ifdef DEBUG
    std::cout<<"---------ShakeA之后--------------"<<std::endl;
    _structure_info_data= DataManager::getInstance().getMDData()->_structure_info_data;
    thrust::copy(_device_data->_d_vx.begin(), _device_data->_d_vx.end(), h_vx.begin());
    thrust::copy(_device_data->_d_px.begin(), _device_data->_d_px.end(), h_px.begin());

    for (int j = 0; j < 10; ++j) {std::cout<<h_vx[j]<<" , ";}
    std::cout<<std::endl;
    for (int j = 0; j < 10; ++j) {std::cout<<h_px[j]<<" , ";}
    std::cout<<std::endl;
#endif

  _force_controller->Execute();
#ifdef DEBUG
    std::cout<<"---------计算力之后--------------"<<std::endl;
    _structure_info_data= DataManager::getInstance().getMDData()->_structure_info_data;
    thrust::copy(_device_data->_d_vx.begin(), _device_data->_d_vx.end(), h_vx.begin());
    thrust::copy(_device_data->_d_px.begin(), _device_data->_d_px.end(), h_px.begin());

    for (int j = 0; j < 10; ++j) {std::cout<<h_vx[j]<<" , ";}
    std::cout<<std::endl;
    for (int j = 0; j < 10; ++j) {std::cout<<h_px[j]<<" , ";}
    std::cout<<std::endl;
#endif
  if ("LANGEVIN"==DataManager::getInstance().getConfigData()->Get<std::string>("temp_ctrl_type", "execution"))
  {
	  _temperature_controller->Update();
  }

  _velocity_controller->Update();
#ifdef DEBUG
    std::cout<<"---------第二次更新速度--------------"<<std::endl;
    _structure_info_data= DataManager::getInstance().getMDData()->_structure_info_data;
    thrust::copy(_device_data->_d_vx.begin(), _device_data->_d_vx.end(), h_vx.begin());
    thrust::copy(_device_data->_d_px.begin(), _device_data->_d_px.end(), h_px.begin());

    for (int j = 0; j < 10; ++j) {std::cout<<h_vx[j]<<" , ";}
    std::cout<<std::endl;
    for (int j = 0; j < 10; ++j) {std::cout<<h_px[j]<<" , ";}
    std::cout<<std::endl;
#endif
  if (use_shake)
  {
    _shake_controller->ShakeB();
  }
#ifdef DEBUG
    std::cout<<"---------ShakeB之后--------------"<<std::endl;
    _structure_info_data= DataManager::getInstance().getMDData()->_structure_info_data;
    thrust::copy(_device_data->_d_vx.begin(), _device_data->_d_vx.end(), h_vx.begin());
    thrust::copy(_device_data->_d_px.begin(), _device_data->_d_px.end(), h_px.begin());

    for (int j = 0; j < 10; ++j) {std::cout<<h_vx[j]<<" , ";}
    std::cout<<std::endl;
    for (int j = 0; j < 10; ++j) {std::cout<<h_px[j]<<" , ";}
    std::cout<<std::endl;
#endif
  _temperature_controller->ComputeTemp();
#ifdef DEBUG
    std::cout<<"---------更新温度--------------"<<std::endl;
    _structure_info_data= DataManager::getInstance().getMDData()->_structure_info_data;
    thrust::copy(_device_data->_d_vx.begin(), _device_data->_d_vx.end(), h_vx.begin());
    thrust::copy(_device_data->_d_px.begin(), _device_data->_d_px.end(), h_px.begin());

    for (int j = 0; j < 10; ++j) {std::cout<<h_vx[j]<<" , ";}
    std::cout<<std::endl;
    for (int j = 0; j < 10; ++j) {std::cout<<h_px[j]<<" , ";}
    std::cout<<std::endl;
#endif
  if ("LANGEVIN" == DataManager::getInstance().getConfigData()->Get<std::string>("temp_ctrl_type", "execution"))
	  return;

  _temperature_controller->Update();

  CHECK_RUNTIME(hipDeviceSynchronize());
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<rbmd::Real> duration = end - start;

  std::cout << "time pre step "<< duration.count() << "秒" << std::endl;
}

void NVTensemble::Postsolve() {}
