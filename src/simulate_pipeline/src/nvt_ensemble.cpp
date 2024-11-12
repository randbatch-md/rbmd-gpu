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
//#define DEBUG
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


    std::shared_ptr<StructureInfoData> _structure_info_data;
    auto _device_data=DataManager::getInstance().getDeviceData();
    _structure_info_data= DataManager::getInstance().getMDData()->_structure_info_data;

    std::vector<rbmd::Real> h_fx(300);
    std::vector<rbmd::Real> h_fy(300);
    std::vector<rbmd::Real> h_fz(300);

    std::vector<rbmd::Real> h_vx(300);
    std::vector<rbmd::Real> h_vy(300);
    std::vector<rbmd::Real> h_vz(300);

    std::vector<rbmd::Real> h_px(300);
    std::vector<rbmd::Real> h_py(300);
    std::vector<rbmd::Real> h_pz(300);

    thrust::copy(_device_data->_d_fx.begin(), _device_data->_d_fx.end(), h_fx.begin());
    thrust::copy(_device_data->_d_fy.begin(), _device_data->_d_fy.end(), h_fy.begin());
    thrust::copy(_device_data->_d_fz.begin(), _device_data->_d_fz.end(), h_fz.begin());

    thrust::copy(_device_data->_d_vx.begin(), _device_data->_d_vx.end(), h_vx.begin());
    thrust::copy(_device_data->_d_vy.begin(), _device_data->_d_vy.end(), h_vy.begin());
    thrust::copy(_device_data->_d_vz.begin(), _device_data->_d_vz.end(), h_vz.begin());

    thrust::copy(_device_data->_d_px.begin(), _device_data->_d_px.end(), h_px.begin());
    thrust::copy(_device_data->_d_py.begin(), _device_data->_d_py.end(), h_py.begin());
    thrust::copy(_device_data->_d_pz.begin(), _device_data->_d_pz.end(), h_pz.begin());

    std::cout<<"导入：力的初始计算结果为："<<std::endl;
    for (int j = 0; j < h_fx.size(); ++j)
    {
        std::cout<<"[ "<<h_fx[j]<<" , "<<h_fy[j]<<" , "<<h_fz[j]<<" ]" <<std::endl;
    }
    std::cout<<std::endl;

    std::cout<<"导入：速度的初始结果为："<<std::endl;
    for (int j = 0; j < h_vx.size(); ++j)
    {
        std::cout<<"[ "<<h_vx[j]<<" , "<<h_vy[j]<<" , "<<h_vz[j]<<" ]" <<std::endl;
    }
    std::cout<<std::endl;

    std::cout<<"导入：位置的初始结果为："<<std::endl;
    for (int j = 0; j < h_px.size(); ++j)
    {
        std::cout<<"[ "<<h_px[j]<<" , "<<h_py[j]<<" , "<<h_pz[j]<<" ]" <<std::endl;
    }
    std::cout<<std::endl;

  _position_controller->Init();
  _velocity_controller->Init();
  _temperature_controller->Init();

  _force_controller->Init();
  _force_controller->Execute();
  _shake_controller->Init();



    thrust::copy(_device_data->_d_fx.begin(), _device_data->_d_fx.end(), h_fx.begin());
    thrust::copy(_device_data->_d_fy.begin(), _device_data->_d_fy.end(), h_fy.begin());
    thrust::copy(_device_data->_d_fz.begin(), _device_data->_d_fz.end(), h_fz.begin());

    thrust::copy(_device_data->_d_vx.begin(), _device_data->_d_vx.end(), h_vx.begin());
    thrust::copy(_device_data->_d_vy.begin(), _device_data->_d_vy.end(), h_vy.begin());
    thrust::copy(_device_data->_d_vz.begin(), _device_data->_d_vz.end(), h_vz.begin());

    thrust::copy(_device_data->_d_px.begin(), _device_data->_d_px.end(), h_px.begin());
    thrust::copy(_device_data->_d_py.begin(), _device_data->_d_py.end(), h_py.begin());
    thrust::copy(_device_data->_d_pz.begin(), _device_data->_d_pz.end(), h_pz.begin());

    //auto h_atoms_id = LinkedCellLocator::GetInstance().GetLinkedCell()->_atom_id_to_idx;
    thrust::host_vector<rbmd::Id> h_atoms_id = _device_data->_d_atoms_id;

    std::cout<<"初始计算：力的初始计算结果为："<<std::endl;
    for (int j = 0; j < h_atoms_id.size(); ++j)
    {
        std::cout<<"[ "<<h_fx[h_atoms_id[j]]<<" , "<<h_fy[h_atoms_id[j]]<<" , "<<h_fz[h_atoms_id[j]]<<" ]" <<std::endl;
    }
    std::cout<<std::endl;

    std::cout<<"初始计算：速度的初始结果为："<<std::endl;
    for (int j = 0; j < h_atoms_id.size(); ++j)
    {
        std::cout<<"[ "<<h_vx[h_atoms_id[j]]<<" , "<<h_vy[h_atoms_id[j]]<<" , "<<h_vz[h_atoms_id[j]]<<" ]" <<std::endl;
    }
    std::cout<<std::endl;

    std::cout<<"初始计算：位置的初始结果为："<<std::endl;
    for (int j = 0; j < h_atoms_id.size(); ++j)
    {
        std::cout<<"[ "<<h_px[h_atoms_id[j]]<<" , "<<h_py[h_atoms_id[j]]<<" , "<<h_pz[h_atoms_id[j]]<<" ]" <<std::endl;
    }
    std::cout<<std::endl;
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
  std::vector<rbmd::Real> h_fx(300);

  thrust::copy(_device_data->_d_vx.begin(), _device_data->_d_vx.end(), h_vx.begin());
  thrust::copy(_device_data->_d_px.begin(), _device_data->_d_px.end(), h_px.begin());

    std::vector<rbmd::Real> h_vy(300);
    std::vector<rbmd::Real> h_vz(300);
    thrust::copy(_device_data->_d_vy.begin(), _device_data->_d_vy.end(), h_vy.begin());
    thrust::copy(_device_data->_d_vz.begin(), _device_data->_d_vz.end(), h_vz.begin());

    for (int j = 0; j < h_vx.size(); ++j)
    {
        std::cout<<"[ "<<h_vx[j]<<" , "<<h_vy[j]<<" , "<<h_vz[j]<<" ]" <<std::endl;
    }

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
