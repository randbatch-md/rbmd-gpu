#include "berendsen_pressure_controller.h"
#include "temperature_controller.h"

#include <thrust/device_ptr.h>

#include <cmath>

#include "device_types.h"
#include "unit_factor.h"
#include "update_pressure_op.h"
#include "common/thermo_stats.hpp"
extern rbmd::Real test_temperature;
extern rbmd::Id test_current_step;
BerendsenPressureController::BerendsenPressureController() {
  std::remove("pressure.txt");
}
BerendsenPressureController::~BerendsenPressureController() {

};

void BerendsenPressureController::Init() {
  auto pressure_array=DataManager::getInstance().getConfigData()->
GetArray<rbmd::Real>("pressure", "execution"); //[1.0,1.0,1.0,10.0]
  _pressure_start = pressure_array[0];
  _pressure_stop = pressure_array[1];
  _pressure_damp = pressure_array[2];
  _bulkmodulus = pressure_array[3];

  _dt =  DataManager::getInstance().getConfigData()->Get<rbmd::Real>(
          "timestep", "execution");//0.001

  auto unit = DataManager::getInstance().getConfigData()->Get
    <std::string>("unit", "init_configuration", "read_data");

  UNIT unit_factor = unit_factor_map[unit];
  switch (unit_factor) {
    case UNIT::LJ:
      _nktv2p = UnitFactor<UNIT::LJ>::_nktv2p;
      _kB = UnitFactor<UNIT::LJ>::_kb;
      break;
    case UNIT::REAL:
      _nktv2p = UnitFactor<UNIT::REAL>::_mvv2e;
      _kB = UnitFactor<UNIT::REAL>::_kb;
      break;
    default:
      break;
  }

  //tdof
  Computedof();
  //
  REAL_DATA(_p_start)[0] = REAL_DATA(_p_start)[1] = REAL_DATA(_p_start)[2]
  = _pressure_start;
   REAL_DATA(_p_stop)[0] = REAL_DATA(_p_stop)[1] = REAL_DATA(_p_stop)[2]
  = _pressure_stop;
   REAL_DATA(_p_damp)[0] = REAL_DATA(_p_damp)[1] = REAL_DATA(_p_damp)[2]
  = _pressure_damp;
  auto temperature_array=DataManager::getInstance().getConfigData()->
  GetArray<rbmd::Real>("temperature", "execution");
  ThermoStats::Instance().AddThermoData("temperature",temperature_array[0]);

  ThermoStats::Instance().AddThermoData("pressure",_pressure_start);
}

void BerendsenPressureController::Update()
{
  ComputePressure();
  Couple();

  //Compute dilation
  auto currentstep = test_current_step;
  auto beginstep = 0;
  auto endstep = DataManager::getInstance().getConfigData()
    ->Get<rbmd::Real>("num_steps", "execution");//10000

  rbmd::Real delta = currentstep - beginstep;
  if (delta != 0.0)
  {
    delta = delta / static_cast<rbmd::Real>(endstep - beginstep);
  }

  for (int i = 0; i < 3; i++)
  {
    auto dt_over_period = _dt / REAL_DATA(_p_damp)[i];
    auto bulkmodulus_inv = 1.0 / _bulkmodulus;
    REAL_DATA(_p_target)[i] = REAL_DATA(_p_start)[i] +
      delta * (REAL_DATA(_p_stop)[i] - REAL_DATA(_p_start)[i]);

    REAL_DATA(_dilation)[i] =
      POW(1.0 - dt_over_period * (REAL_DATA(_p_target)[0] -
        REAL_DATA(_p_current)[i])* bulkmodulus_inv, 1.0 / 3.0);
  }

  //reset box and atoms
  ResetBox();
}

void BerendsenPressureController::ComputePressure()
{
  auto volume = CalculateVolume(*_box);
  auto  inv_volume = 1/volume;

  ComputeVirial();

  //compute pressure
  _pressure = (_tdof * _kB * test_temperature+_device_data->_d_virial[0]
    + _device_data->_d_virial[1] +_device_data->_d_virial[2])
  /3.0 * inv_volume * _nktv2p;

  ThermoStats::Instance().AddThermoData("pressure",_pressure);
  // out
  std::ofstream outfile("pressure.txt", std::ios::app);
  outfile << test_current_step << " " << _pressure << std::endl;
  outfile.close();
}

void BerendsenPressureController::ComputeVirial()
{
  TransformForces(_device_data->_d_virial,_device_data->_d_virial_lj,
    _device_data->_d_virial_specialcoul,_device_data->_d_virial_kspace,
    _device_data->_d_virial_bond,_device_data->_d_virial_angle,
    _device_data->_d_virial_dihedral);
}
void BerendsenPressureController::Computedof()
{
  auto num_atoms = *(_structure_info_data->_num_atoms);
  auto extra_dof = 3; //dimension =3
  _tdof = 3 * num_atoms - extra_dof;
  bool use_shake = false;
  if (use_shake)
  {
    _tdof = _tdof - num_atoms;
  }
}

void BerendsenPressureController::Couple()
{
  REAL_DATA(_p_current)[0] = REAL_DATA(_p_current)[1] =
  REAL_DATA(_p_current)[2] = _pressure;
}

void BerendsenPressureController::ResetBox()
{
  //  convert lamda coords
  X2Lamda();

  // change range and box

  rbmd::Real oldlo, oldhi, ctr;
  bool pbc[3] = {1, 1, 1};
  for (int i = 0; i < 3; i++)
  {
    oldlo = _box->_coord_min[i];
    oldhi = _box->_coord_max[i];
    ctr = 0.5 * (oldlo + oldhi);
    _box->_coord_min[i] = (oldlo - ctr) * REAL_DATA(_dilation)[i] + ctr;
    _box->_coord_max[i] = (oldhi - ctr) * REAL_DATA(_dilation)[i] + ctr;
  }
  _box->Setup(_box->_type, _box->_coord_min, _box->_coord_max, pbc);

  // CHECK_RUNTIME(
  //     MEMCPY(_box., h_box, sizeof(Box), H2D));

  // std::cout << "range.Min=" <<  _box->_coord_min[0] << ",range.Max="
  // << _box->_coord_max[0] << std::endl;

  // convert real coords
  Lamda2X();

}

void BerendsenPressureController::X2Lamda(){
  auto num_atoms = *(_structure_info_data->_num_atoms);

  op::X2LamdaOp<device::DEVICE_GPU>()(
    *_box,num_atoms,
    thrust::raw_pointer_cast(_device_data->_d_px.data()),
    thrust::raw_pointer_cast(_device_data->_d_py.data()),
    thrust::raw_pointer_cast(_device_data->_d_pz.data()));
}

void BerendsenPressureController::Lamda2X(){
  auto num_atoms = *(_structure_info_data->_num_atoms);

  op::Lamda2XOp<device::DEVICE_GPU>()(
    *_box,num_atoms,
    thrust::raw_pointer_cast(_device_data->_d_px.data()),
    thrust::raw_pointer_cast(_device_data->_d_py.data()),
    thrust::raw_pointer_cast(_device_data->_d_pz.data()));
}


