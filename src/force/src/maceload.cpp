#include <torch/torch.h>
#include <torch/script.h>
#include <algorithm>
#include <iostream>
#include <array>
#include "../include/maceload.h"
#include <torch/jit.h>
#include <chrono>
#include <thrust/transform.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include "../../data_manager/include/data_manager.h"
#include "../../data_manager/include/model/device_data.h"
#include "../../data_manager/include/model/box.h"
#include "../../data_manager/include/scheduler/memory_scheduler.h"
#include "../../data_manager/include/model/structure_info_data.h"
#include "neighbor_list/include/neighbor_list_builder/full_neighbor_list_builder.h"
#include "../../src/neighbor_list/src/op/full_neighbor_list_op.h"
//maceload macetest;
maceload::maceload(std::shared_ptr<DeviceData>& device_data,
          const std::shared_ptr<StructureData>& structure_data,
          const std::shared_ptr<StructureInfoData>& structure_info_data)
    : _device_data(device_data),  // 初始化引用
      _structure_data(structure_data),  // 初始化常量引用
      _structure_info_data(structure_info_data)   {
  //_device_data = std::make_shared<DeviceData>();
  _neighbor_list_builder = std::make_shared<FullNeighborListBuilder>();
  _box = std::make_shared<Box>();

};
void maceload::Init()
{
  torch::jit::getProfilingMode() = false;
  torch::jit::getExecutorMode() = true; 
  auto& num_atoms = *(_structure_info_data->_num_atoms);
  n_nodes = num_atoms;
  std::cout << n_nodes<<" atoms " << std::endl;
  //输出原子数
  vflag_global = false;
  std::string path_ = "./mace_agnesi_small.model-lammps.pt";
  //可优化
  //if (device_name != "cuda" && device_name != "dcu")
  //{
  //  device = c10::Device(torch::kCPU);
  //}
  try {
        model = torch::jit::load(path_, device);
        std::cout << "Loading MACE model from " << path_ << "\" ..."<< std::endl;
    }
  catch (const c10::Error& e) {
        std::cerr << "Error load mace model\n";
    }

  mace_r_max = model.attr("r_max").toTensor().item<double>();
    //double num_inter = model.attr("num_interactions").toTensor().item<double>();
    //std::cout << num_inter <<std::endl;
  auto mace_atom_table = model.attr("atomic_numbers").toTensor();
  n_node_feats = mace_atom_table.numel();
  std::string transsymbol;
  int transnum;
  float transmass;
  for (int a_n = 0; a_n < n_node_feats; ++a_n) {
        transnum = mace_atom_table[a_n].item<int>() -1;
        transmass = mass_periodic_table[transnum];
        mace_feats_mass_table.push_back(transmass);
        transsymbol = periodic_table[transnum];
        mace_feats_table.push_back(transsymbol);
        //std::cout<< transsymbol << std::endl;
        //std::cout << transmass << std::endl;
  }    
  energy = energy.to(device);
  batch = torch::zeros({ n_nodes }, torch::dtype(torch::kInt64).device(device));
  forces = torch::empty({ n_nodes, 3 }, torch_float_dtype.device(device));
  ptr[1] = n_nodes;
  ptr = ptr.to(device);
  weight[0] = 1.0;
  weight = weight.to(device);
  mask = torch::ones({ n_nodes }, torch::dtype(torch::kBool));
  node_attrs = torch::zeros({ n_nodes, n_node_feats }, torch_float_dtype.device(device));
  const rbmd::Real* d_mass = thrust::raw_pointer_cast(_device_data->_d_mass.data());;
  for (int n_a_i = 0; n_a_i < n_nodes; ++n_a_i)
  {
    
    float am = d_mass[n_a_i];
    float dvalue;
    int dvalueflag = 0;
    for (int n_f_i = 0; n_f_i < n_node_feats; ++n_f_i)
    {
      dvalue = am - mace_feats_mass_table[n_f_i];
      if (dvalue < 0.05 && dvalue > -0.05)
      {
        node_attrs[n_a_i][n_f_i] = 1.0;
        dvalueflag++;
      }
    }
    if (dvalueflag != 1)
    {
      std::cout << "No." << n_a_i
                << " atom have an inaccurate atomic mass.\n Please check or enter atomic symbol."
                << std::endl;
    }
  }
}
void maceload::Execute()
{
  loadcell();
  loadpositions();
  loadedges();
  forward();
}

void maceload::loadatoms(std::vector<std::string> symbollist) {
    /* no change between two timestep */
    batch = torch::zeros({ n_nodes }, torch::dtype(torch::kInt64).device(device));
    forces = torch::empty({ n_nodes, 3 }, torch_float_dtype.device(device));

    ptr[1] = n_nodes;
    ptr = ptr.to(device);
    weight[0] = 1.0;
    weight = weight.to(device);
    node_attrs = torch::zeros({ n_nodes, n_node_feats }, torch_float_dtype.device(device));

    for (int n_a_i = 0; n_a_i < n_nodes; ++n_a_i) {
        //std::string as = symbollist[n_a_i];
        std::string as = symbollist[n_a_i];
        //大小写问题
        try {
            auto iter = std::find(mace_feats_table.begin(), mace_feats_table.end(), as);
            node_attrs[n_a_i][std::distance(mace_feats_table.begin(), iter)] = 1.0;
        }
        catch (const c10::Error& e) {
            std::cerr << "Unsupported element symbol!\n";
        }
    }
    mask = torch::ones({ n_nodes }, torch::dtype(torch::kBool));
}
void maceload::loadmass(thrust::device_vector<rbmd::Real> masslist)
{
  /* no change between two timestep */
  batch = torch::zeros({ n_nodes }, torch::dtype(torch::kInt64).device(device));
  forces = torch::empty({ n_nodes, 3 }, torch_float_dtype.device(device));
  ptr[1] = n_nodes;
  ptr = ptr.to(device);
  weight[0] = 1.0;
  weight = weight.to(device);
  mask = torch::ones({ n_nodes }, torch::dtype(torch::kBool));

  node_attrs = torch::zeros({ n_nodes, n_node_feats }, torch_float_dtype.device(device));
  for (int n_a_i = 0; n_a_i < n_nodes; ++n_a_i)
  {
    float am = masslist[n_a_i];
    float dvalue;
    int dvalueflag = 0;
    for (int n_f_i = 0; n_f_i < n_node_feats; ++n_f_i)
    {
      dvalue = am - mace_feats_mass_table[n_f_i];
      if (dvalue < 0.05 && dvalue > -0.05)
      {
        node_attrs[n_a_i][n_f_i] = 1.0;
        dvalueflag++;
      }
    }
    if (dvalueflag != 1)
    {
      std::cout << "No." << n_a_i
                << " atom have an inaccurate atomic mass.\n Please check or enter atomic symbol."
                << std::endl;
    }
  }
}
void maceload::loadcell()
{
    #if USE_DOUBLE  
    cell[0][0] = _box->_length[0];
    cell[1][1] = _box->_length[1];
    cell[2][2] = _box->_length[2];
    cell[2][1] = _box->_length[3];
    cell[2][0] = _box->_length[4];
    cell[1][0] = _box->_length[5];
    #else
    cell[0][0] = static_cast<double>(_box->_length[0]);
    cell[1][1] = static_cast<double>(_box->_length[1]);
    cell[2][2] = static_cast<double>(_box->_length[2]);
    cell[2][1] = static_cast<double>(_box->_length[3]);
    cell[2][0] = static_cast<double>(_box->_length[4]);
    cell[1][0] = static_cast<double>(_box->_length[5]);
    #endif
    /*
    h[0] = xprd;
    h[1] = yprd;
    h[2] = zprd;
    h_inv[0] = 1.0/h[0];
    h_inv[1] = 1.0/h[1];
    h_inv[2] = 1.0/h[2];
    h[3] = yz;
    h[4] = xz;
    h[5] = xy;
    cell[2][1] = h[3];
    cell[2][0] = h[4];
    cell[1][0] = h[5];
    h_inv[3] = -h[3] / (h[1]*h[2]);
    h_inv[4] = (h[3]*h[5] - h[1]*h[4]) / (h[0]*h[1]*h[2]);
    h_inv[5] = -h[5] / (h[0]*h[1]);*/
    cell = cell.to(device);
}
void maceload::loadpositions()
{
    const rbmd::Real* d_px = thrust::raw_pointer_cast(_device_data->_d_px.data());
    const rbmd::Real* d_py = thrust::raw_pointer_cast(_device_data->_d_py.data());
    const rbmd::Real* d_pz = thrust::raw_pointer_cast(_device_data->_d_pz.data());
    #if USE_DOUBLE
    positions = torch::empty({ n_nodes, 3 }, torch::kFloat64).to(torch::kCUDA);
    thrust::device_ptr<rbmd::Real> ptr_tensor = thrust::device_pointer_cast(positions.data_ptr<rbmd::Real>());
    thrust::copy(d_px, d_px + n_nodes, ptr_tensor);                  
    thrust::copy(d_py, d_py + n_nodes, ptr_tensor + n_nodes);           
    thrust::copy(d_pz, d_pz + n_nodes, ptr_tensor + n_nodes*2); 
    #else
    positions = torch::empty({ n_nodes, 3 }, torch::kFloat32).to(torch::kCUDA);
    thrust::device_ptr<rbmd::Real> ptr_tensor = thrust::device_pointer_cast(positions.data_ptr<rbmd::Real>());
    thrust::copy(d_px, d_px + n_nodes, ptr_tensor);                  
    thrust::copy(d_py, d_py + n_nodes, ptr_tensor + n_nodes);          
    thrust::copy(d_pz, d_pz + n_nodes, ptr_tensor + n_nodes*2);  
    positions = positions.to(torch::kF64);
    #endif
}
void loadedges_index_8(rbmd::Id* per_atom_cell_id,
            rbmd::Id* in_atom_list_start_index,
            rbmd::Id* in_atom_list_end_index, rbmd::Real cutoff_2,
            rbmd::Id total_atom_num, rbmd::Real* px, rbmd::Real* py,
            rbmd::Real* pz, rbmd::Id* max_neighbor_num,
            rbmd::Id* neighbor_start, rbmd::Id* neighbor_end,
            rbmd::Id* neighbors,rbmd::Id* neighbors_atoms,
            rbmd::Real* unit_shiftx,rbmd::Real* unit_shifty,
            rbmd::Real* unit_shiftz,rbmd::Real* shiftx,
            rbmd::Real* shifty,rbmd::Real* shiftz,
            Box* d_box, rbmd::Id* should_realloc,
            rbmd::Id* neighbor_cell, rbmd::Id neighbor_cell_num) {



}
void maceload::loadedges()
{
    //std::cout << "edges_start" << std::endl;
    //auto start_time1 = std::chrono::high_resolution_clock::now();
    //std::cout << index << std::endl;
  _list = _neighbor_list_builder->Build();
  long index = _device_data->_d_neighbors_atoms.size();
  rbmd::Id* d_neighbor = thrust::raw_pointer_cast(_neighbor_list->_d_neighbors.data());
  rbmd::Id* d_neighbor_atoms = thrust::raw_pointer_cast(_device_data->_d_neighbors_atoms.data());
  const rbmd::Real* d_unit_shiftx = thrust::raw_pointer_cast(_device_data->_d_unit_shiftx.data());
  const rbmd::Real* d_unit_shifty = thrust::raw_pointer_cast(_device_data->_d_unit_shifty.data());
  const rbmd::Real* d_unit_shiftz = thrust::raw_pointer_cast(_device_data->_d_unit_shiftz.data());
  const rbmd::Real* d_shiftx = thrust::raw_pointer_cast(_device_data->_d_shiftx.data());
  const rbmd::Real* d_shifty = thrust::raw_pointer_cast(_device_data->_d_shifty.data());
  const rbmd::Real* d_shiftz = thrust::raw_pointer_cast(_device_data->_d_shiftz.data());

    #if USE_DOUBLE
  unit_shifts =
torch::cat({ torch::from_blob(d_unit_shiftx, {index , 1}, optsfloat64),
             torch::from_blob(d_unit_shifty, {  index ,1}, optsfloat64),
             torch::from_blob(d_unit_shiftz, {index , 1}, optsfloat64)},
           1);
  shifts =
    torch::cat({ torch::from_blob(d_shiftx, {index , 1}, optsfloat64),
                 torch::from_blob(d_shifty, {  index ,1}, optsfloat64),
                 torch::from_blob(d_shiftz, {index , 1}, optsfloat64)},
               1);
    #else
  unit_shifts = torch::empty({ index, 3 }, torch::kFloat32).to(torch::kCUDA);
  thrust::device_ptr<rbmd::Real> ustr_tensor = thrust::device_pointer_cast(unit_shifts.data_ptr<rbmd::Real>());
  thrust::copy(d_unit_shiftx, d_unit_shiftx + index, ustr_tensor);
  thrust::copy(d_unit_shifty, d_unit_shifty + index, ustr_tensor + index);
  thrust::copy(d_unit_shiftz, d_unit_shiftz + index, ustr_tensor + index*2);
  unit_shifts = unit_shifts.to(torch::kF64);
  shifts = torch::empty({ index, 3 }, torch::kFloat32).to(torch::kCUDA);
  thrust::device_ptr<rbmd::Real> str_tensor = thrust::device_pointer_cast(shifts.data_ptr<rbmd::Real>());
  thrust::copy(d_shiftx, d_shiftx + index, str_tensor);
  thrust::copy(d_shifty, d_shifty + index, str_tensor + index);
  thrust::copy(d_shiftz, d_shiftz + index, str_tensor + index*2);
  shifts = shifts.to(torch::kF64);
#endif

    #if USE_64BIT_IDS
  edge_index = torch::cat({ torch::from_blob(d_neighbor, { 1, index }, optsint64),
             torch::from_blob(d_neighbor_atoms, { 1, index }, optsint64) },
           0);

    #else

    edge_index =
        torch::cat({ torch::from_blob(d_neighbor, { 1, index }, optsint32),
                     torch::from_blob(d_neighbor_atoms, { 1, index }, optsint32) },
                   0)
          .to(torch::kInt64);
    #endif


}

double maceload::energyout(c10::impl::GenericDict ouput)
{
  energy = ouput.at("total_energy_local").toTensor();
  eng_vdwl = energy.item<double>();
  return eng_vdwl;
}

void maceload::forward() {
    //std::cout << "forward_start"<< std::endl;
    //auto start_time = std::chrono::high_resolution_clock::now();
    c10::Dict<std::string, torch::Tensor> input;
    try
    {
      edge_index = edge_index.to(device);
      //positions = positions.to(device);
      shifts = shifts.to(device);
      unit_shifts = unit_shifts.to(device);
    }
    catch (const c10::Error& e)
    {
      std::cerr << "Error in to device\n";
    }
    //std::cout << edge_index.sizes() << std::endl;
    input.insert("batch", batch);
    input.insert("cell", cell);
    input.insert("edge_index", edge_index);
    input.insert("energy", energy);
    input.insert("forces", forces);
    input.insert("node_attrs", node_attrs);
    input.insert("positions", positions);
    input.insert("ptr", ptr);
    input.insert("shifts", shifts);
    input.insert("unit_shifts", unit_shifts);
    input.insert("weight", weight);
    //std::cout << "数据载入完成" << std::endl;
    model.eval();
    c10::impl::GenericDict output = model.forward({ input, mask.to(device), bool(vflag_global) }).toGenericDict();
    //auto end_time = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    //std::cout << "Time elapsed: " << duration.count() << " milliseconds" << std::endl;
    //std::cout << "模型前向完成" << std::endl;
    #if USE_DOUBLE
    energy = output.at("total_energy_local").toTensor();//gpu double
    double eng_vdwl = energy.item<double>();
    std::cout << eng_vdwl<< std::endl;
    forces = output.at("forces").toTensor();//gpu float
    double* forces_data_ptr = forces.data_ptr<double>();
    thrust::copy(forces_data_ptr, forces_data_ptr + n_nodes, _device_data->_d_fx.begin());     // X方向
    thrust::copy(forces_data_ptr + n_nodes, forces_data_ptr + 2 * n_nodes, _device_data->_d_fy.begin());  // Y方向
    thrust::copy(forces_data_ptr + 2 * n_nodes, forces_data_ptr + 3 * n_nodes, _device_data->_d_fz.begin()); // Z方向
    if (vflag_global) 
    {
        vir = output.at("virials").toTensor();
        _device_data->_d_virial_xx[0] += vir[0][0][0].item<double>();
        _device_data->_d_virial_yy[0] += vir[0][1][1].item<double>();
        _device_data->_d_virial_zz[0] += vir[0][2][2].item<double>();
        _device_data->_d_virial_xy[0] += 0.5 * (vir[0][1][0].item<double>() + vir[0][0][1].item<double>());
        _device_data->_d_virial_xz[0] += 0.5 * (vir[0][2][0].item<double>() + vir[0][0][2].item<double>());
        _device_data->_d_virial_yz[0] += 0.5 * (vir[0][2][1].item<double>() + vir[0][1][2].item<double>());
    }
    #else
    energy = output.at("total_energy_local").toTensor().to(torch::kFloat);//gpu double
    float eng_vdwl = energy.item<float>();
    std::cout << eng_vdwl<< std::endl;
    forces = output.at("forces").toTensor().to(torch::kFloat);//gpu float
    float* forces_data_ptr = forces.data_ptr<float>();
    thrust::copy(forces_data_ptr, forces_data_ptr + n_nodes, _device_data->_d_fx.begin());     // X方向
    thrust::copy(forces_data_ptr + n_nodes, forces_data_ptr + 2 * n_nodes, _device_data->_d_fy.begin());  // Y方向
    thrust::copy(forces_data_ptr + 2 * n_nodes, forces_data_ptr + 3 * n_nodes, _device_data->_d_fz.begin()); // Z方向
    if (vflag_global) 
    {
        vir = output.at("virials").toTensor().to(torch::kFloat);
        _device_data->_d_virial_xx[0] += vir[0][0][0].item<float>();
        _device_data->_d_virial_yy[0] += vir[0][1][1].item<float>();
        _device_data->_d_virial_zz[0] += vir[0][2][2].item<float>();
        _device_data->_d_virial_xy[0] += 0.5 * (vir[0][1][0].item<float>() + vir[0][0][1].item<float>());
        _device_data->_d_virial_xz[0] += 0.5 * (vir[0][2][0].item<float>() + vir[0][0][2].item<float>());
        _device_data->_d_virial_yz[0] += 0.5 * (vir[0][2][1].item<float>() + vir[0][1][2].item<float>());
    }
    #endif
}

torch::Tensor maceload::virialout(c10::impl::GenericDict ouput) {
    torch::Tensor vir = ouput.at("virials").toTensor().cpu();
    torch::Tensor virial = torch::zeros({6}, torch_float_dtype);
    virial[0] += vir[0][0][0].item<double>();
    virial[1] += vir[0][1][1].item<double>();
    virial[2] += vir[0][2][2].item<double>();
    virial[3] += 0.5 * (vir[0][1][0].item<double>() + vir[0][0][1].item<double>());
    virial[4] += 0.5 * (vir[0][2][0].item<double>() + vir[0][0][2].item<double>());
    virial[5] += 0.5 * (vir[0][2][1].item<double>() + vir[0][1][2].item<double>());
    return virial;
}
torch::Tensor maceload::atomenergy(c10::impl::GenericDict ouput) {
    torch::Tensor node_energy = ouput.at("node_energy").toTensor();
    return node_energy;
}