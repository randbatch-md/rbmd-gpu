#pragma once
#include <thrust/host_vector.h>

#include <memory>

#include "data_manager.h"
#include "model/device_data.h"
#include "model/structure_info_data.h"

class Force {
 public:
  Force() {
    this->_structure_info_data =
        DataManager::getInstance().getMDData()->_structure_info_data;
    this->_device_data = DataManager::getInstance().getDeviceData();
    this->_box = DataManager::getInstance().getMDData()->_box;
  };
  virtual ~Force() = default;

  // virtual void Update()=0;
  virtual void Init() {};
  virtual void Execute() = 0;
  virtual void EvaluatePotentialenergy(){};
  void ReduceVirial(
    rbmd::Id num_atoms,
    const  thrust::device_vector<rbmd::Real>& d_flat_virial_atom,
    thrust::device_vector<rbmd::Real>& d_virial)
  {
    thrust::host_vector<rbmd::Real> h_flat_virial_atom(d_flat_virial_atom);

    std::vector<rbmd::Real> virial(6, 0.0);
    for(int atom = 0; atom < num_atoms; ++atom){
      for(int j = 0; j < 6; ++j){
        virial[j] += h_flat_virial_atom[j * num_atoms + atom];
      }
    }

    //H2D
    thrust::copy(virial.begin(),
    virial.end(), d_virial.begin());
  }

 protected:
  std::shared_ptr<StructureInfoData> _structure_info_data;
  std::shared_ptr<DeviceData> _device_data;

  std::shared_ptr<Box> _box;
};