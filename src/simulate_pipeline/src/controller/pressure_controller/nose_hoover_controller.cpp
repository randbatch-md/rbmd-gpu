#include "nose_hoover_controller.h"

// #include <thrust/device_ptr.h>
//
// #include <cmath>
//
#include "common/thermo_stats.hpp"
#include "default_position_controller.h"
#include "default_velocity_controller.h"
// #include "device_types.h"
#include "../simulate_pipeline/src/controller/group_controller/group_controller.h"
// #include "simulate.h"
//
#include "unit_factor.h"
#include "update_pressure_op.h"
#include "update_temperature_op.h"
extern int test_current_step;

NoseHooverController::NoseHooverController() {
  std::remove("temperature.txt");
  CHECK_RUNTIME(MALLOC(&_d_temp_contrib, sizeof(rbmd::Real)));

  _position_controller = std::make_shared<DefaultPositionController>();
  _velocity_controller = std::make_shared<DefaultVelocityController>();

  _group_controller = std::make_shared<GroupController>();
}
NoseHooverController::~NoseHooverController() {
  CHECK_RUNTIME(FREE(_d_temp_contrib));
};

void NoseHooverController::Init() {

  _velocity_controller->Init();
  _position_controller->Init();

  _ensemble_type =DataManager::getInstance().getConfigData()->
  Get<std::string>("ensemble", "execution");
  if("NPT" == _ensemble_type) {
    _pressure_flag = true;
  }

  //read target temperature_array
  auto temperature_array=DataManager::getInstance().getConfigData()->
    GetArray<rbmd::Real>("temperature", "execution"); //[1.0,1.0,1.0]
  _t_start = temperature_array[0];
  _t_stop = temperature_array[1];
  _t_damp = temperature_array[2];

  //read target pressure_array
  if(_pressure_flag) {
    auto pressure_array=DataManager::getInstance().getConfigData()->
      GetArray<rbmd::Real>("pressure", "execution"); //[1.0,1.0,1.0,10.0]
    _pressure_start = pressure_array[0];
    _pressure_stop = pressure_array[1];
    _pressure_damp = pressure_array[2];
    _bulkmodulus = pressure_array[3];
  }

  //read timestep
  _dt =  DataManager::getInstance().getConfigData()->Get<rbmd::Real>(
          "timestep", "execution");//0.001

  auto unit = DataManager::getInstance().getConfigData()->Get
    <std::string>("unit", "init_configuration", "read_data");

  UNIT unit_factor = unit_factor_map[unit];
  switch (unit_factor) {
    case UNIT::LJ:
      _nktv2p = UnitFactor<UNIT::LJ>::_nktv2p;
      _mvv2e = UnitFactor<UNIT::LJ>::_mvv2e;
      _kB = UnitFactor<UNIT::LJ>::_kb;
      _fmt2v = UnitFactor<UNIT::LJ>::_fmt2v;
      break;
    case UNIT::METAL:
      _nktv2p = UnitFactor<UNIT::METAL>::_nktv2p;
      _mvv2e = UnitFactor<UNIT::METAL>::_mvv2e;
      _kB = UnitFactor<UNIT::METAL>::_kb;
      _fmt2v = UnitFactor<UNIT::METAL>::_fmt2v;
      break;
    case UNIT::REAL:
      _nktv2p = UnitFactor<UNIT::REAL>::_nktv2p;
      _mvv2e = UnitFactor<UNIT::REAL>::_mvv2e;
      _kB = UnitFactor<UNIT::REAL>::_kb;
      _fmt2v = UnitFactor<UNIT::REAL>::_fmt2v;
      break;
    default:
      break;
  }

  //compute tdof
  Computedof();

  //pressure init
  REAL_DATA(_p_start)[0] = REAL_DATA(_p_start)[1] = REAL_DATA(_p_start)[2]
  = _pressure_start;
  REAL_DATA(_p_stop)[0]  = REAL_DATA(_p_stop)[1]  = REAL_DATA(_p_stop)[2]
  = _pressure_stop;
  REAL_DATA(_p_damp)[0]  = REAL_DATA(_p_damp)[1]  = REAL_DATA(_p_damp)[2]
  = _pressure_damp;
  REAL_DATA(_p_freq)[0] = REAL_DATA(_p_freq)[1] = REAL_DATA(_p_freq)[2]
  = 1 / _pressure_damp;

  _t_freq = 1 / _t_damp;

  _p_freq_max = 0.0;
  _p_freq_max = MAX(REAL_DATA(_p_freq)[0], REAL_DATA(_p_freq)[1]);
  _p_freq_max = MAX(_p_freq_max, REAL_DATA(_p_freq)[2]);

  _p_flag.resize(6);
  _p_flag[0] = _p_flag[1] = _p_flag[2] = 1;
  _pdim = _p_flag[0] + _p_flag[1] + _p_flag[2];

  //defult values
  _eta_mass_flag = 1;
  _nc_tchain = _nc_pchain = 1;
  _mtchain = _mpchain = 3;  //defult value =3
  _mtk_flag = 1;

  // set timesteps
  _dtf = _dt * _fmt2v;
  _dthalf = 0.5 * _dt;
  _dt4 = 0.25 * _dt;
  _dt8 = 0.125 * _dt;
  _dto = _dthalf;

  _drag = 0.0;
  _tdrag_factor = 1.0 - (_dt * _t_freq * _drag / _nc_tchain);
  _pdrag_factor = 1.0 - (_dt * _p_freq_max * _drag / _nc_pchain);

  // Nose-Hoover thermostat  init
  _eta.resize(_mtchain);
  _eta_dot.resize(_mtchain+1);
  _eta_dot[_mtchain] = 0;
  _eta_dotdot.resize(_mtchain);
  _eta_mass.resize(_mtchain);

  for (int ich = 0; ich < _mtchain; ich++)
  {
    _eta[ich] = _eta_dot[ich] = _eta_dotdot[ich] = 0.0;
  }

  //Nose-Hoover barostat init
  _omega.resize(6);
  _omega_dot.resize(6);
  _omega_mass.resize(6);

  _omega[0] = _omega[1] = _omega[2] = 0.0;
  _omega_dot[0] = _omega_dot[1] = _omega_dot[2] = 0.0;
  _omega_mass[0] = _omega_mass[1] = _omega_mass[2] = 0.0;

  _omega[3] = _omega[4] = _omega[5] = 0.0;
  _omega_dot[3] = _omega_dot[4] = _omega_dot[5] = 0.0;
  _omega_mass[3] = _omega_mass[4] = _omega_mass[5] = 0.0;

  _etap.resize(_mtchain);
  _etap_dot.resize(_mpchain + 1);
  _etap_dot[_mpchain] = 0.0;
  _etap_dotdot.resize(_mtchain);
  _etap_mass.resize(_mtchain);
  for (int ich = 0; ich < _mpchain; ich++)
  {
    _etap[ich] = _etap_dot[ich] = _etap_dotdot[ich] = 0.0;
  }

  //读取group并初始化GroupController
  const auto& config = DataManager::getInstance().getConfigData();
  if (config->PathExists({"execution","group"})) {
    _group_name = DataManager::getInstance().getConfigData()->Get<std::string>(
    "group", "execution");
    if (_group_name.empty()) {
      _group_name = "all";
    }
    _group_controller->Init();
  }


  //com_bias
  if (config->PathExists({"execution","com_bias"}))
  {
    auto com_bias = DataManager::getInstance().getConfigData()->Get<std::string>(
     "com_bias", "execution");
    if ("yes" == com_bias ) {
      _com_bias = true;
    }
  }

  //Nose-Hoover parameters init
  SetUp();

  ThermoStats::Instance().AddThermoData("temperature",_t_start);
  ThermoStats::Instance().AddThermoData("pressure",_pressure_start);
}

void NoseHooverController::Update()
{
}

void NoseHooverController::ComputeTemperature(){
  rbmd::Id num_atoms = *(_structure_info_data->_num_atoms);
  CHECK_RUNTIME(MEMSET(_d_temp_contrib, 0, sizeof(rbmd::Real)));

  if (_com_bias) {
    //  计算质心速度 (vbias)
   _group_controller->ComputeVCM(_group_name, _vbias);
    op::ComputeTemperatureCOMOp<device::DEVICE_GPU>()(num_atoms, _mvv2e,_vbias,
        thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
        thrust::raw_pointer_cast(_device_data->_d_mass.data()),
        thrust::raw_pointer_cast(_device_data->_d_vx.data()),
        thrust::raw_pointer_cast(_device_data->_d_vy.data()),
        thrust::raw_pointer_cast(_device_data->_d_vz.data()), _d_temp_contrib);
  }
  else {
    op::ComputeTemperatureOp<device::DEVICE_GPU>()(num_atoms, _mvv2e,
    thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
    thrust::raw_pointer_cast(_device_data->_d_mass.data()),
    thrust::raw_pointer_cast(_device_data->_d_vx.data()),
    thrust::raw_pointer_cast(_device_data->_d_vy.data()),
    thrust::raw_pointer_cast(_device_data->_d_vz.data()), _d_temp_contrib);
  }


  CHECK_RUNTIME(MEMCPY(&_temp_sum, _d_temp_contrib, sizeof(rbmd::Real), D2H));

  _temperature = 0.5 * _temp_sum / (_tdof * _kB / 2.0);

  if (std::isnan(_temperature)) {
    Logger::Instance().error( "\033[31mFATAL ERROR: The temperature of the MD simulation is NaN"
                             ". Please check the initial model and the force field parameters. "
    "is invalid.\033[0m");
    exit(EXIT_FAILURE); //
  }
}

void NoseHooverController::ComputeVirial()
{
  //need to modify....
  TransformForces(_device_data->_d_virial,_device_data->_d_virial_lj,
    _device_data->_d_virial_kspace,_device_data->_d_virial_bond,
    _device_data->_d_virial_angle,_device_data->_d_virial_dihedral);
}

void NoseHooverController::ComputePressure()
{
  auto volume = CalculateVolume(*_box);
  auto  inv_volume = 1/volume;

  ComputeVirial();

  //compute pressure
  _pressure = (_tdof * _kB * _temperature+_device_data->_d_virial[0]
    + _device_data->_d_virial[1] +_device_data->_d_virial[2])
  /3.0 * inv_volume * _nktv2p;

}

void NoseHooverController::Couple()
{
  REAL_DATA(_p_current)[0] = REAL_DATA(_p_current)[1] =
    REAL_DATA(_p_current)[2] = _pressure;
}

void NoseHooverController::SetUp()
{
  rbmd::Id num_atoms = *(_structure_info_data->_num_atoms);

  ComputeTemperature();       // current temperature
  ComputeTempTarget();  //target temperature and ke

  if (_pressure_flag)
  {
    ComputePressTarget();  //target pressure
    ComputePressure();   //current pressure
    Couple();
  }

  // masses and initial forces on thermostat variables
  _eta_mass[0] = _tdof * _kB * _t_target / (_t_freq * _t_freq);
  for (int ich = 1; ich < _mtchain; ich++)
  {
    _eta_mass[ich] = _kB * _t_target / (_t_freq * _t_freq);
  }

  for (int ich = 1; ich < _mtchain; ich++)
  {
    _eta_dotdot[ich] = (_eta_mass[ich - 1] * _eta_dot[ich - 1] *
      _eta_dot[ich - 1] - _kB * _t_target) /_eta_mass[ich];
  }

  if (_pressure_flag)
  {
    // masses and initial forces on barostat variables
    rbmd::Real kt = _kB * _t_target;
    rbmd::Real nkt = (num_atoms + 1) * kt;
    for (int i = 0; i < 3; i++)
    {
      if (_p_flag[i])
      {
        _omega_mass[i] = nkt / (REAL_DATA(_p_freq)[i] * REAL_DATA(_p_freq)[i]);
      }
    }
    // masses and initial forces on barostat thermostat variables
    if (_mpchain)
    {
      _etap_mass[0] = _kB* _t_target / (_p_freq_max * _p_freq_max);
      for (int ich = 1; ich < _mpchain; ich++)
      {
        _etap_mass[ich] = _kB * _t_target / (_p_freq_max * _p_freq_max);
      }

      for (int ich = 1; ich < _mpchain; ich++)
      {
        _etap_dotdot[ich] = (_etap_mass[ich - 1] * _etap_dot[ich - 1] * _etap_dot[ich - 1] -
                            _kB * _t_target) /_etap_mass[ich];
      }
    }
  }

}

void NoseHooverController::ComputeTempTarget()
{
    if (_t_stop == _t_start) //Thermostatic simulation
    {
        _t_target = _t_stop= _t_start;
    }
    else                    //anisothermal simulation
    {
        auto currentstep = test_current_step;
        auto beginstep = 0;
        auto endstep = DataManager::getInstance().getConfigData()->
          Get<rbmd::Real>("num_steps", "execution");

        rbmd::Real delta = currentstep - beginstep;

        if (delta != 0.0)
        {
            delta = delta / static_cast<rbmd::Real>(endstep - beginstep);
        }

        _t_target = _t_start + delta * (_t_stop - _t_start);
    }
    //
    _ke_target = _tdof * _kB * _t_target;
}

void NoseHooverController::ComputePressTarget()
{
    _p_hydro = 0.0;
    for (int i = 0; i < 3; i++)
    {
        if (REAL_DATA(_p_stop)[i] == REAL_DATA(_p_start)[i]) //Isobaric simulation
        {
            REAL_DATA(_p_target)[i] = REAL_DATA(_p_stop)[i] = REAL_DATA(_p_start)[i];
        }
        else
        {
          auto currentstep = test_current_step;
          auto beginstep = 0;
          auto endstep = DataManager::getInstance().getConfigData()->
            Get<rbmd::Real>("num_steps", "execution");

            rbmd::Real delta = currentstep - beginstep;
            if (delta != 0.0)
            {
                delta = delta / static_cast<rbmd::Real>(endstep - beginstep);
            }
            for (int i = 0; i < 3; i++)
            {
                REAL_DATA(_p_target)[i] = REAL_DATA(_p_start)[i] + delta *
                  (REAL_DATA(_p_stop)[i] - REAL_DATA(_p_start)[i]);
            }
        }

        //
        _p_hydro += REAL_DATA(_p_target)[i];
        if (_pdim > 0)
        {
            _p_hydro /= _pdim;
        }
    }
 //TRICLINIC TODO:
  // if deviatoric, recompute sigma each time p_target changes
}

void NoseHooverController::NHOmegaDot()
{
  rbmd::Id num_atoms = *(_structure_info_data->_num_atoms);

  //
  rbmd::Real f_omega;
  rbmd::Real volume = CalculateVolume(*_box);


  //mtk_term1
  _mtk_term1 = 0.0;
  if (_mtk_flag)
  {
    _mtk_term1 = _tdof * _kB * _temperature;
    _mtk_term1 /= _pdim * num_atoms;
  }
  for (int i = 0; i < 3; i++)
  {
    if (_p_flag[i])
    {
      f_omega = (REAL_DATA(_p_current)[i] - _p_hydro) * volume / (_omega_mass[i] * _nktv2p)
      + _mtk_term1 / _omega_mass[i];
      _omega_dot[i] += f_omega * _dthalf;
      _omega_dot[i] *= _pdrag_factor;
    }
  }

  //mtk_term2
  _mtk_term2 = 0.0;
  if (_mtk_flag)
  {
    for (int i = 0; i < 3; i++)
    {
      if (_p_flag[i])
      {
        _mtk_term2 += _omega_dot[i];
      }
    }
    if (_pdim > 0)
    {
      _mtk_term2 /= _pdim * num_atoms;
    }
  }
}

void NoseHooverController::NH_V_Press()
{
  Real3 factor;
  REAL_DATA(factor)[0] = EXP(-_dt4 * (_omega_dot[0] + _mtk_term2));
  REAL_DATA(factor)[1] = EXP(-_dt4 * (_omega_dot[1] + _mtk_term2));
  REAL_DATA(factor)[2] = EXP(-_dt4 * (_omega_dot[2] + _mtk_term2));

  // perform half-step barostat scaling of velocities
  if (_com_bias) {
    //1: remove-bias
    op::RemoveBiasOp<device::DEVICE_GPU>()(
        *(_structure_info_data->_num_atoms),_vbias,
        thrust::raw_pointer_cast(_device_data->_d_vx.data()),
        thrust::raw_pointer_cast(_device_data->_d_vy.data()),
        thrust::raw_pointer_cast(_device_data->_d_vz.data()));
    //2 : scale
    op::UpdataVelocityRescalePressureOp<device::DEVICE_GPU>()(
        *(_structure_info_data->_num_atoms), factor,
             thrust::raw_pointer_cast(_device_data->_d_vx.data()),
             thrust::raw_pointer_cast(_device_data->_d_vy.data()),
             thrust::raw_pointer_cast(_device_data->_d_vz.data()));

    op::UpdataVelocityRescalePressureOp<device::DEVICE_GPU>()(
    *(_structure_info_data->_num_atoms), factor,
         thrust::raw_pointer_cast(_device_data->_d_vx.data()),
         thrust::raw_pointer_cast(_device_data->_d_vy.data()),
         thrust::raw_pointer_cast(_device_data->_d_vz.data()));
    //3: restore-bias
    op::RestoreBiasOp<device::DEVICE_GPU>()(
    *(_structure_info_data->_num_atoms),_vbias,
    thrust::raw_pointer_cast(_device_data->_d_vx.data()),
    thrust::raw_pointer_cast(_device_data->_d_vy.data()),
    thrust::raw_pointer_cast(_device_data->_d_vz.data()));
  }
  else {
    op::UpdataVelocityRescalePressureOp<device::DEVICE_GPU>()(
            *(_structure_info_data->_num_atoms), factor,
                 thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                 thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                 thrust::raw_pointer_cast(_device_data->_d_vz.data()));

    op::UpdataVelocityRescalePressureOp<device::DEVICE_GPU>()(
               *(_structure_info_data->_num_atoms), factor,
                    thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                    thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                    thrust::raw_pointer_cast(_device_data->_d_vz.data()));
  }

}

void NoseHooverController::InitialIntegrate()
{
  if (_pressure_flag)
  {
    // update eta_press_dot
    NHCPressIntegrate();
  }

  // update eta_dot
  ComputeTempTarget();  // target temperature
  NHCTempIntegrate();  //perform half-step update of chain thermostat variables

  // need to recompute pressure to account for change in KE
  if (_pressure_flag)
  {
    //ComputeTempe();     //current temperature
    ComputePressure();  //current pressure
    Couple();

     //
    ComputePressTarget(); // target pressure
    NHOmegaDot();
    NH_V_Press();
  }

  _velocity_controller->Update();

  if (_pressure_flag)
  {
    ResetBox();    // reset box in the first half-step
  }

  _position_controller->Update();

  if (_pressure_flag)
  {
    ResetBox(); // Reset the box in the second half-step
  }
}

void NoseHooverController::FinalIntegrate()
{

  _velocity_controller->Update();

  if (_pressure_flag)
  {
    NH_V_Press();
  }

  // need to compute new temperature and pressure after velocities rescaled
  ComputeTemperature(); // current temperature

  if (_pressure_flag)
  {
    ComputePressure(); // current pressure
    Couple();
    NHOmegaDot();
  }

  // update eta_dot
  // update eta_press_dot
  NHCTempIntegrate();

  if (_pressure_flag)
  {
    NHCPressIntegrate();
  }

  //
  ThermoStats::Instance().AddThermoData("temperature",_temperature);
  ThermoStats::Instance().AddThermoData("pressure",_pressure);

  //out
  auto interval = DataManager::getInstance().getConfigData()->Get<rbmd::Id>(
"interval", "outputs", "thermo_out");
  std::ofstream outfile("temperature.txt", std::ios::app);
  if (outfile.tellp() == 0) {
    outfile << "step temperature pressure" << std::endl;
  }
  if (test_current_step % interval == 0) {
    outfile << test_current_step << " " << _temperature  << " "<< _pressure
  << std::endl;
  }
  outfile.close();
}

void NoseHooverController::NHCTempIntegrate()
{
  rbmd::Real expfac;
  rbmd::Real ke_current = _tdof * _kB * _temperature;

  // Update masses, to preserve initial freq, if flag set
  if (_eta_mass_flag)
  {
    _eta_mass[0] = _tdof * _kB * _t_target / (_t_freq * _t_freq);
    for (int ich = 1; ich < _mtchain; ich++)
    {
      _eta_mass[ich] = _kB * _t_target / (_t_freq * _t_freq);
    }
  }

  if (_eta_mass[0] > 0.0)
  {
    _eta_dotdot[0] = (ke_current - _ke_target) / _eta_mass[0]; //链的加速度
  }else{
    _eta_dotdot[0] = 0.0;
  }

  //
  rbmd::Real ncfac = 1.0 / _nc_tchain;
  for (int iloop = 0; iloop < _nc_tchain; iloop++)
  {
    for (int ich = _mtchain - 1; ich > 0; ich--) //This must be done starting from the last link and proceeding forward to the first link.
                                                //This ensures that the updates of the upstream links (those closer to the particles)
                                                //can accurately reflect the influence of the downstream links.

    {
      expfac = EXP(-ncfac * _dt8 * _eta_dot[ich + 1]);  //The velocity is updated exponentially using exp(-dt/8),
                                             //and the ith link segment will be influenced by the subsequent link segment 𝜂(ih + 1).
      _eta_dot[ich] *= expfac;
      _eta_dot[ich] += _eta_dotdot[ich] * ncfac * _dt4;  // Based on the current link acceleration eta_dotdot,
                                                      //  the speed of the chain is updated. The time step is dt/4.
      _eta_dot[ich] *= _tdrag_factor;                  //
      _eta_dot[ich] *= expfac;                        // Once again, use exp(-dt/8) for exponential update to complete another half-step update.
                                                     //This step ensures the complete update of momentum, making the evolution conform to the time-reversal symmetry.
    }                                                //The application of the exponential scaling factor twice is
                                                          //to simulate the Hamiltonian equations of the Nose-Hoover chain system.

    expfac = EXP(-ncfac * _dt8 * _eta_dot[1]);
    _eta_dot[0] *= expfac;
    _eta_dot[0] += _eta_dotdot[0] * ncfac * _dt4;
    _eta_dot[0] *= _tdrag_factor;
    _eta_dot[0] *= expfac;

    _factor_eta = EXP(-ncfac * _dthalf * _eta_dot[0]);

    //nh_v_temp();  Update of dt/2
    if (_com_bias) {
      //1: remove-bias
      op::RemoveBiasOp<device::DEVICE_GPU>()(
          *(_structure_info_data->_num_atoms),_vbias,
          thrust::raw_pointer_cast(_device_data->_d_vx.data()),
          thrust::raw_pointer_cast(_device_data->_d_vy.data()),
          thrust::raw_pointer_cast(_device_data->_d_vz.data()));
      //2 : scale
      op::UpdataVelocityRescaleOp<device::DEVICE_GPU>()(
          *(_structure_info_data->_num_atoms), _factor_eta,
          thrust::raw_pointer_cast(_device_data->_d_vx.data()),
          thrust::raw_pointer_cast(_device_data->_d_vy.data()),
          thrust::raw_pointer_cast(_device_data->_d_vz.data()));
      //3: restore-bias
      op::RestoreBiasOp<device::DEVICE_GPU>()(
        *(_structure_info_data->_num_atoms),_vbias,
        thrust::raw_pointer_cast(_device_data->_d_vx.data()),
        thrust::raw_pointer_cast(_device_data->_d_vy.data()),
        thrust::raw_pointer_cast(_device_data->_d_vz.data()));
    }
    else {
      op::UpdataVelocityRescaleOp<device::DEVICE_GPU>()(
        *(_structure_info_data->_num_atoms), _factor_eta,
        thrust::raw_pointer_cast(_device_data->_d_vx.data()),
        thrust::raw_pointer_cast(_device_data->_d_vy.data()),
        thrust::raw_pointer_cast(_device_data->_d_vz.data()));
    }
    // rescale temperature due to velocity scaling
     _temperature = _temperature * _factor_eta * _factor_eta;

    //updata  ke
    ke_current = _tdof * _kB * _temperature;

    if (_eta_mass[0] > 0.0) {
      _eta_dotdot[0] = (ke_current - _ke_target) / _eta_mass[0];
    }else {
      _eta_dotdot[0] = 0.0;
    }

    for (int ich = 0; ich < _mtchain; ich++) {
      _eta[ich] += ncfac * _dthalf * _eta_dot[ich];   //The position of the chain
    }

    _eta_dot[0] *= expfac;
    _eta_dot[0] += _eta_dotdot[0] * ncfac * _dt4;
    _eta_dot[0] *= expfac;

    for (int ich = 1; ich < _mtchain; ich++) //Forward loop (from beginning to end): Used for updating acceleration
    {
      expfac = EXP(-ncfac * _dt8 * _eta_dot[ich + 1]);
      _eta_dot[ich] *= expfac;
      _eta_dotdot[ich] =(_eta_mass[ich - 1] * _eta_dot[ich - 1] *
        _eta_dot[ich - 1] - _kB * _t_target) / _eta_mass[ich];
      _eta_dot[ich] += _eta_dotdot[ich] * ncfac * _dt4;
      _eta_dot[ich] *= expfac;
    }
  }
}

void NoseHooverController::NHCPressIntegrate()
{
  rbmd::Id pdof;
  rbmd::Real expfac, factor_etap, ke_current;
  rbmd::Real kt = _kB * _t_target;
  rbmd::Real lkt_press;

  ke_current = 0.0;
  pdof = 0;
  for (int i = 0; i < 3; i++)
  {
    if (_p_flag[i])
    {
      ke_current += _omega_mass[i] * _omega_dot[i] * _omega_dot[i];
      pdof++;
    }
  }

  //
  lkt_press = kt;  //iso
  _etap_dotdot[0] = (ke_current - lkt_press) / _etap_mass[0]; // dotdot :  the pressure bath chain

  //
  rbmd::Real ncfac = 1.0 / _nc_pchain;
  for (int iloop = 0; iloop < _nc_pchain; iloop++)
  {

    for (int ich = _mpchain - 1; ich > 0; ich--) //counterpropagation
    {
      expfac = EXP(-ncfac * _dt8 * _etap_dot[ich + 1]);
      _etap_dot[ich] *= expfac;
      _etap_dot[ich] += _etap_dotdot[ich] * ncfac * _dt4;
      _etap_dot[ich] *= _pdrag_factor;
      _etap_dot[ich] *= expfac;
    }

    expfac = EXP(-ncfac * _dt8 * _etap_dot[1]);
    _etap_dot[0] *= expfac;
    _etap_dot[0] += _etap_dotdot[0] * ncfac * _dt4;
    _etap_dot[0] *= _pdrag_factor;
    _etap_dot[0] *= expfac;              //The speed of the pressure bath chain

    for (int ich = 0; ich < _mpchain; ich++)
    {
      _etap[ich] += ncfac * _dthalf * _etap_dot[ich];
    }


    factor_etap = EXP(-ncfac * _dthalf * _etap_dot[0]);
    for (int i = 0; i < 3; i++)
    {
      if (_p_flag[i])
      {
        _omega_dot[i] *= factor_etap; //ISO
      }
    }

    ke_current = 0.0;
    for (int i = 0; i < 3; i++)
    {
      if (_p_flag[i])
      {
        ke_current += _omega_mass[i] * _omega_dot[i] * _omega_dot[i];
      }
    }

    _etap_dotdot[0] = (ke_current - lkt_press) / _etap_mass[0];

    _etap_dot[0] *= expfac;
    _etap_dot[0] += _etap_dotdot[0] * ncfac * _dt4;
    _etap_dot[0] *= expfac;

    for (int ich = 1; ich < _mpchain; ich++) //Forward loop (from beginning to end): Used for updating acceleration
    {
      expfac = EXP(-ncfac * _dt8 * _etap_dot[ich + 1]);
      _etap_dot[ich] *= expfac;
      _etap_dotdot[ich] = (_etap_mass[ich - 1] * _etap_dot[ich - 1] *
        _etap_dot[ich - 1] -_kB * _t_target) / _etap_mass[ich];
      _etap_dot[ich] += _etap_dotdot[ich] * ncfac * _dt4;
      _etap_dot[ich] *= expfac;
    }
  }
}

void NoseHooverController::ResetBox()
{
  rbmd::Real oldlo, oldhi;
  rbmd::Real expfac;

  rbmd::Id num_atoms = *(_structure_info_data->_num_atoms);

  //  convert lamda coords
  X2Lamda();

  //double dto2 = dto / 2.0;
  //double dto4 = dto / 4.0;
  //double dto8 = dto / 8.0;

  if (_p_flag[0])
  {
    oldlo = _box->_coord_min[0];
    oldhi = _box->_coord_max[0];
    expfac = EXP(_dto * _omega_dot[0]);
    _box->_coord_min[0] = (oldlo - _box->_median_point[0]) * expfac
    + _box->_median_point[0];
    _box->_coord_max[0] = (oldhi - _box->_median_point[0]) * expfac
    + _box->_median_point[0];
  }

  if (_p_flag[1])
  {
    oldlo = _box->_coord_min[1];
    oldhi = _box->_coord_max[1];
    expfac = EXP(_dto * _omega_dot[1]);
    _box->_coord_min[1] = (oldlo - _box->_median_point[1]) * expfac
    + _box->_median_point[1];
    _box->_coord_max[1]  = (oldhi - _box->_median_point[1]) * expfac
    + _box->_median_point[1];
    //if (scalexy)
    //{
    //  h[5] *= expfac;
    //}
  }

  if (_p_flag[2])
  {
    oldlo = _box->_coord_min[2];
    oldhi = _box->_coord_max[2];
    expfac = EXP(_dto * _omega_dot[2]);
    _box->_coord_min[2] = (oldlo - _box->_median_point[2]) * expfac
    + _box->_median_point[2];
    _box->_coord_max[2]= (oldhi - _box->_median_point[2]) * expfac
    + _box->_median_point[2];
    //if (scalexz)
    //{
    //  h[4] *= expfac;
    //}
    //if (scaleyz)
    //{
    //  h[3] *= expfac;
    //}
  }

  bool pbc[3] = {1, 1, 1};
  _box->Setup(_box->_type, _box->_coord_min, _box->_coord_max, pbc);

  // convert real coords
  Lamda2X();

}

void NoseHooverController::X2Lamda(){
  auto num_atoms = *(_structure_info_data->_num_atoms);

  op::X2LamdaOp<device::DEVICE_GPU>()(
    *_box,num_atoms,
    thrust::raw_pointer_cast(_device_data->_d_px.data()),
    thrust::raw_pointer_cast(_device_data->_d_py.data()),
    thrust::raw_pointer_cast(_device_data->_d_pz.data()));
}

void NoseHooverController::Lamda2X(){
  auto num_atoms = *(_structure_info_data->_num_atoms);

  op::Lamda2XOp<device::DEVICE_GPU>()(
    *_box,num_atoms,
    thrust::raw_pointer_cast(_device_data->_d_px.data()),
    thrust::raw_pointer_cast(_device_data->_d_py.data()),
    thrust::raw_pointer_cast(_device_data->_d_pz.data()));
}

void NoseHooverController::Computedof()
{
  rbmd::Id num_atoms = *(_structure_info_data->_num_atoms);
  auto extra_dof = 3; //dimension =3
  _tdof = 3 * num_atoms - extra_dof;
  bool shake = DataManager::getInstance().getConfigData()->GetJudge<bool>
  ( "fix_shake", "hyper_parameters", "extend");
  if (shake == true)
  {
    _tdof = _tdof - num_atoms;
  }
}

void NoseHooverController::InitialBM() {
    // ---------------------------------------------------------------------------------
    //  步骤 1: 恒温器/恒压器演化 & 第一次速度/位置缩放 (Operator A: for dt/2)
    // ---------------------------------------------------------------------------------
    // 演化恒温器链 (thermostat chain) for dt/2
    ComputeTempTarget();
    NHCTempIntegrate(); // 该函数应包含半步积分逻辑 (内部使用 _dthalf)

    // 如果是NPT，演化恒压器链 (barostat chain) for dt/2
    if (_pressure_flag)
    {
      //ComputeTempe();     //current temperature
      ComputePressure();  //current pressure
      Couple();

      //
      ComputePressTarget(); // target pressure
      NHOmegaDot();
      NH_V_Press();
    }

    if (_pressure_flag)
    {
      ResetBox();    // reset box in the first half-step
    }

    // ---------------------------------------------------------------------------------
    //  步骤 2: 核心 Beeman 积分 (Operator B: for dt)
    // ---------------------------------------------------------------------------------
    // 2a. Beeman 位置更新
    _position_controller->Updatebm();

    if (_pressure_flag)
    {
      ResetBox(); // Reset the box in the second half-step
    }
}

void NoseHooverController::FinalBM() {
  // 2c. Beeman 速度更新 (校正步)
  // 使用新力 F(t+dt), 中间力 F(t), 和上一步的力 F(t-dt)
  _velocity_controller->Updatebm();

  if (_pressure_flag)
  {
    NH_V_Press();
  }

  // need to compute new temperature and pressure after velocities rescaled
  ComputeTemperature(); // current temperature

  if (_pressure_flag)
  {
    ComputePressure(); // current pressure
    Couple();
    NHOmegaDot();
  }

  // update eta_dot
  // update eta_press_dot
  NHCTempIntegrate();

  if (_pressure_flag)
  {
    NHCPressIntegrate();
  }

  // ---------------------------------------------------------------------------------
  //  步骤 3: 输出热力学统计信息
  // ---------------------------------------------------------------------------------
  ThermoStats::Instance().AddThermoData("temperature", _temperature);
  if (_pressure_flag) {
    ThermoStats::Instance().AddThermoData("pressure", _pressure);
  }
}
