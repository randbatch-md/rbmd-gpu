#include <algorithm>
#include <iostream>
#include <iomanip>
#include <array>
#include <c10/util/Logging.h>
#include <torch/jit.h>
#include <thrust/binary_search.h>
#include <thrust/scan.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/execution_policy.h>
#include <chrono>
#include <thrust/transform.h>
#include <thrust/device_vector.h>
#include <thrust/gather.h>
#include <thrust/copy.h>
#include <thrust/sequence.h>
#include <thrust/tuple.h>
#include "../../data_manager/include/data_manager.h"
#include "../../data_manager/include/model/device_data.h"
#include "../../data_manager/include/model/box.h"
#include "../../data_manager/include/scheduler/memory_scheduler.h"
#include "../../data_manager/include/model/structure_info_data.h"
#include "neighbor_list/include/neighbor_list_builder/mace_neighbor_list_builder.h"
#include "../../src/neighbor_list/src/op/mace_neighbor_list_op.h"
#include "neighbor_list/include/linked_cell/linked_cell_locator.h"
#include "../include/maceload.h"

#include <spdlog/sinks/basic_file_sink.h>
#include <torch/torch.h>
#include <fstream>
#include <thrust/device_ptr.h>


#undef REDUCE
#include <torch/script.h>
#define REDUCE hipcub::DeviceReduce::Sum

extern int test_current_step;
//maceload macetest;
//maceload::maceload(std::shared_ptr<DeviceData>& device_data,          const std::shared_ptr<StructureData>& structure_data,          const std::shared_ptr<StructureInfoData>& structure_info_data)
maceload::maceload()
{
  //_device_data = std::make_shared<DeviceData>();
  //_mace_device_data = std::make_shared<DeviceData>();
  //_structure_data = std::make_shared<StructureData>();
  //_mace_structure_info_data = std::make_shared<StructureInfoData>();
  _neighbor_list_builder = std::make_shared<MACENeighborListBuilder>();
  _mace_device_data= DataManager::getInstance().getDeviceData();
  _linked_cell = LinkedCellLocator::GetInstance().GetLinkedCell();
  //_box = std::make_shared<Box>();
  std::remove("thermo.txt");
};
void maceload::Init()
{
  torch::jit::getProfilingMode() = false;
  torch::jit::getExecutorMode() = true;
  auto& num_atoms = *(_structure_info_data->_num_atoms);
  n_nodes = num_atoms;//不能正常读入
  #if USE_MACE
    std::cout << "ENABLED" << std::endl;
  #endif
  std::cout << n_nodes<<" atoms " << std::endl;
  //输出原子数
  vflag_global = false;
  std::string path_ = DataManager::getInstance().getConfigData()->Get<std::string>("model_file", "hyper_parameters", "force_field");
  //std::string path_ = "./mace_agnesi_small.model-lammps.pt";
  //可优化
  //if (device_name != "cuda" && device_name != "dcu")
  //{
  //  device = c10::Device(torch::kCPU);
  //}
  if (!torch::cuda::is_available()) {
    std::cout << "CUDA unavailable, setting device type to torch::kCPU." << std::endl;
    device = c10::Device(torch::kCPU);
  } else {
    std::cout << "CUDA found, setting device type to torch::kCUDA." << std::endl;
    device = c10::Device(torch::kCUDA,0);
  }
  try {
        model = torch::jit::load(path_, device);
        std::cout << "Loading MACE model from " << path_ << "\" ..."<< std::endl;
    }
  catch (const c10::Error& e) {
        std::cerr << "Error load mace model\n";
    }

  r_max = model.attr("r_max").toTensor();
  mace_r_max = r_max.item<double>();
  std::cout<<"MACE cutoff :"<<mace_r_max<<" Angstrom"<<std::endl;
  _linked_cell->_cutoff = mace_r_max+0.0001;
  //_linked_cell->_cutoff= 7.0;
  //double num_inter = model.attr("num_interactions").toTensor().item<double>();
  //初始化不同model的mace_table
  //std::cout<<"1"<<std::endl;
  auto mace_atom_table = model.attr("atomic_numbers").toTensor();
  n_node_feats = mace_atom_table.numel();
  std::string transsymbol;
  int transnum;
  float transmass;
  //std::cout<<"2"<<std::endl;
  for (int a_n = 0; a_n < n_node_feats; ++a_n) {
        transnum = mace_atom_table[a_n].item<int>() -1;
        transmass = mass_periodic_table[transnum];
        mace_feats_mass_table.push_back(transmass);
        transsymbol = periodic_table[transnum];
        mace_feats_table.push_back(transsymbol);
  }
  auto model_type =r_max.dtype();
  if(model_type==torch::kFloat64) {
    torch_float_dtype = torch::dtype(torch::kF64);
  }
  else if(model_type==torch::kFloat32) {
    torch_float_dtype = torch::dtype(torch::kF32);
  }
  //std::cout<<"3"<<std::endl;
  cell = torch::zeros({ 3,3 }, torch_float_dtype);
  weight = torch::empty({ 1 }, torch_float_dtype);
  energy = torch::empty({ 1 }, torch_float_dtype);
  energy = energy.to(device);
  batch = torch::zeros({ n_nodes }, torch_int_dtype.device(device));
  forces = torch::empty({ n_nodes, 3 }, torch_float_dtype.device(device));
  ptr = torch::zeros({ 2 }, torch_int_dtype);
  ptr[1] = n_nodes;
  ptr = ptr.to(device);
  weight[0] = 1.0;
  weight = weight.to(device);
  mask = torch::ones({ n_nodes }, torch::dtype(torch::kBool)).to(device);

  const rbmd::Real* d_mass = thrust::raw_pointer_cast(_mace_device_data->_d_mass.data());
  //std::cout<<"1 2"<<std::endl;
  size_t dmasssize = _mace_device_data->_d_mass.size();
 /* for (size_t i = 0; i < _mace_device_data->_d_atoms_type.size(); ++i) {
    std::cout << "Element " << i << ": " << _mace_device_data->_d_atoms_type[i] << std::endl;
  }*/
  //std::cout<<"2"<< dmasssize<<" 2"<<std::endl;
  std::vector<rbmd::Real> h_mass(dmasssize);
  cudaMemcpy(h_mass.data(), d_mass, dmasssize * sizeof(rbmd::Real), cudaMemcpyDeviceToHost);
  int dvalueflag;
  float amass;
  float dvalue;
  int atoms_flag;
  for (int n_a_i = 0; n_a_i < dmasssize; ++n_a_i)
  {
    //std::cout<<"st"<<std::endl;
    amass = h_mass[n_a_i];
    dvalueflag = 0;
    //std::cout<<"2 2"<<std::endl;
    for (int n_f_i = 0; n_f_i < n_node_feats; ++n_f_i)
    {
      dvalue = amass - mace_feats_mass_table[n_f_i];
      if (dvalue < 0.05 && dvalue > -0.05)
      {
        atoms_flag=n_f_i;
        //node_attrs[n_a_i][n_f_i] = 1.0;
        dvalueflag++;
      }
    }
    if (dvalueflag != 1)
    {
      std::cout<<dvalueflag<<std::endl;
      std::cout << "No." << n_a_i
                << " atom have an inaccurate atomic mass.\n Please check or enter atomic symbol."
                << std::endl;
    }
    else {
      atoms_table.push_back(atoms_flag);
    }
  }
  //std::cout<<"1 3"<<std::endl;
  d_atoms_table = torch::tensor(atoms_table, torch::kLong).to(device).contiguous();
  options = torch::TensorOptions().dtype(torch::kInt).device(device);
  node_attrs = torch::zeros({n_nodes, n_node_feats},
                              torch_float_dtype).to(device);
  //load_atomid();
  model.eval();
  //loadedges();
  std::cout <<"init over" << std::endl;
}
void maceload::Execute()
{
  //std::cout <<"edgesb"<<std::endl;
  auto start = std::chrono::high_resolution_clock::now();
  loadedges();
  auto edges = std::chrono::high_resolution_clock::now();
  loadcell();
  loadpositions();
  load_atomid();
  //std::cout <<"id"<<std::endl;
  forward();
  auto end = std::chrono::high_resolution_clock::now();
  auto neighbor = std::chrono::duration_cast<std::chrono::milliseconds>(edges- start);
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - edges);
  //std::cout<<"edges time: "<<neighbor.count()<<" ms"<<std::endl;
  std::cout<<"forward time: "<<duration.count()<<" micros"<<std::endl;
  //TimingStatistics::Instance().record("Short-Range",duration_rbl_force.count());

}
void maceload::load_atomid() {
  //auto start_time = std::chrono::high_resolution_clock::now();
  //std::cout<<"id start"<<std::endl;
  node_attrs=torch::zeros_like(node_attrs);
  //auto start_time2 = std::chrono::high_resolution_clock::now();
  int* atoms_type_ptr = thrust::raw_pointer_cast(_device_data->_d_atoms_type.data());
  //auto start_time3 = std::chrono::high_resolution_clock::now();
  torch::Tensor atom_indices = torch::from_blob(atoms_type_ptr,
                                     {n_nodes},
                                     options);
  //auto start_time4 = std::chrono::high_resolution_clock::now();
  torch::Tensor mapped_types = d_atoms_table.index({atom_indices});

  //auto start_time5 = std::chrono::high_resolution_clock::now();
  node_attrs.scatter_(1, mapped_types.unsqueeze(1), 1.0f);
  //node_attrs.index_put_({torch::arange(n_nodes, torch::kInt64).to(device),  // 行索引
    // mapped_types},                                    // 列索引
    //torch::ones(1, torch_float_dtype).to(device)       // 填充值
    //);

  /*auto end_time = std::chrono::high_resolution_clock::now();
  auto d1 = std::chrono::duration_cast<std::chrono::microseconds>(start_time2 - start_time);
  auto d2 = std::chrono::duration_cast<std::chrono::microseconds>(start_time3 - start_time2);
  auto d3 = std::chrono::duration_cast<std::chrono::microseconds>(start_time4 - start_time3);
  auto d4 = std::chrono::duration_cast<std::chrono::microseconds>(start_time5- start_time4);
  auto d5 = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time5);
  std::cout << "d1: " << d1.count() << " microseconds"
  << " d2: " << d2.count() << " microseconds"
  << "d3: "<< d3.count() << " microseconds"
  << "d4: " << d4.count() << " microseconds"
  << "d5: "<< d5.count() << " microseconds"<< std::endl;*/

}
void maceload::loadatoms(std::vector<std::string> symbollist) {
    /* no change between two timestep */
  auto model_type =r_max.dtype();
  if(model_type==torch::kFloat64) {
    torch_float_dtype = torch::dtype(torch::kF64);
  }
  else if(model_type==torch::kFloat32) {
    torch_float_dtype = torch::dtype(torch::kF32);
  }
  cell = torch::zeros({ 3,3 }, torch_float_dtype);
  weight = torch::empty({ 1 }, torch_float_dtype);
  energy = torch::empty({ 1 }, torch_float_dtype);
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
  auto model_type =r_max.dtype();
  if(model_type==torch::kFloat64) {
    torch_float_dtype = torch::dtype(torch::kF64);
  }
  else if(model_type==torch::kFloat32) {
    torch_float_dtype = torch::dtype(torch::kF32);
  }
  cell = torch::zeros({ 3,3 }, torch_float_dtype);
  weight = torch::empty({ 1 }, torch_float_dtype);
  energy = torch::empty({ 1 }, torch_float_dtype);
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
    auto model_type =r_max.dtype();
    #if USE_DOUBLE
    if(model_type==torch::kFloat64) {
      cell[0][0] = _box->_length[0];
      cell[1][1] = _box->_length[1];
      cell[2][2] = _box->_length[2];
      cell[2][1] = _box->_length[3];
      cell[2][0] = _box->_length[4];
      cell[1][0] = _box->_length[5];
    }
    else if (model_type==torch::kFloat32) {
      cell[0][0] = static_cast<float>(_box->_length[0]);
      cell[1][1] = static_cast<float>(_box->_length[1]);
      cell[2][2] = static_cast<float>(_box->_length[2]);
      cell[2][1] = static_cast<float>(_box->_length[3]);
      cell[2][0] = static_cast<float>(_box->_length[4]);
      cell[1][0] = static_cast<float>(_box->_length[5]);
    }
    #else
    if(model_type==torch::kFloat64) {
      cell[0][0] = static_cast<double>(_box->_length[0]);
      cell[1][1] = static_cast<double>(_box->_length[1]);
      cell[2][2] = static_cast<double>(_box->_length[2]);
      cell[2][1] = static_cast<double>(_box->_length[3]);
      cell[2][0] = static_cast<double>(_box->_length[4]);
      cell[1][0] = static_cast<double>(_box->_length[5]);
    }
    else if (model_type==torch::kFloat32) {
      cell[0][0] = _box->_length[0];
      cell[1][1] = _box->_length[1];
      cell[2][2] = _box->_length[2];
      cell[2][1] = _box->_length[3];
      cell[2][0] = _box->_length[4];
      cell[1][0] = _box->_length[5];
    }
    /*cell[0][0] = static_cast<double>(_device_data->_d_box->_length[0]);
    cell[1][1] = static_cast<double>(_mace_device_data->_d_box->_length[1]);
    cell[2][2] = static_cast<double>(_mace_device_data->_d_box->_length[2]);
    cell[2][1] = static_cast<double>(_mace_device_data->_d_box->_length[3]);
    cell[2][0] = static_cast<double>(_mace_device_data->_d_box->_length[4]);
    cell[1][0] = static_cast<double>(_mace_device_data->_d_box->_length[5]);*/
    //std::cout <<cell <<std::endl;
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
  //thrust::device_vector<rbmd::Real> d_px = _mace_device_data->_d_px;
  //thrust::device_vector<rbmd::Real> d_py = _mace_device_data->_d_py;
  //thrust::device_vector<rbmd::Real> d_pz = _mace_device_data->_d_pz;
  /*positions = torch::empty({ n_nodes,3}, torch::kFloat64).to(torch::kCUDA);
  for (int i=0; i<n_nodes;++i) {
    positions[i][0]=double(_mace_device_data->_d_px[i]);
    positions[i][1]=double(_mace_device_data->_d_py[i]);
    positions[i][2]=double(_mace_device_data->_d_pz[i]);
  }*/
    /*rbmd::Real d_px = _mace_device_data->_d_px;
    rbmd::Real* d_py = thrust::raw_pointer_cast(_mace_device_data->_d_py.data());
    rbmd::Real* d_pz = thrust::raw_pointer_cast(_mace_device_data->_d_pz.data());*/
  //使用 指针
    //auto start_time = std::chrono::high_resolution_clock::now();
    auto model_type =r_max.dtype();
    const rbmd::Real* d_px = thrust::raw_pointer_cast(_device_data->_d_px.data());
    const rbmd::Real* d_py = thrust::raw_pointer_cast(_device_data->_d_py.data());
    const rbmd::Real* d_pz = thrust::raw_pointer_cast(_device_data->_d_pz.data());
    #if USE_DOUBLE
    positions = torch::empty({ 3 ,n_nodes}, torch::kFloat64).to(torch::kCUDA);
    thrust::device_ptr<rbmd::Real> ptr_tensor = thrust::device_pointer_cast(positions.data_ptr<rbmd::Real>());
    thrust::copy(d_px, d_px + n_nodes, ptr_tensor);
    thrust::copy(d_py, d_py + n_nodes, ptr_tensor + n_nodes);
    thrust::copy(d_pz, d_pz + n_nodes, ptr_tensor + n_nodes*2);
  if(model_type==torch::kFloat64) {
    positions =positions.transpose(0, 1).contiguous();
  }
  else if(model_type==torch::kFloat32) {
    positions =positions.transpose(0, 1).to(torch::kF32).contiguous();
  }
    #else
  //positions = torch::empty({ n_nodes, 3 }, torch::kFloat32).to(torch::kCUDA);
  // 使用直接的内存指针填充tensor
  // 每一列依次填充
  //auto end_time_zzytensor = std::chrono::high_resolution_clock::now();
  //positions.index({torch::indexing::Slice(), 0}) = torch::from_blob(d_px, {n_nodes}, torch::kCUDA);
  //positions.index({torch::indexing::Slice(), 1}) = torch::from_blob(d_py, {n_nodes}, torch::kCUDA);
  //positions.index({torch::indexing::Slice(), 2}) = torch::from_blob(d_pz, {n_nodes}, torch::kCUDA);
  //torch自带，耗时长
  //下方为利用指针直接运行的，但需要转置
    positions = torch::empty({ 3 ,n_nodes}, torch::kFloat32).to(torch::kCUDA);
    //auto end_time_zzytensor = std::chrono::high_resolution_clock::now();
    thrust::device_ptr<rbmd::Real> ptr_tensor = thrust::device_pointer_cast(positions.data_ptr<rbmd::Real>());
    thrust::copy(d_px, d_px + n_nodes, ptr_tensor);                  
    thrust::copy(d_py, d_py + n_nodes, ptr_tensor + n_nodes);          
    thrust::copy(d_pz, d_pz + n_nodes, ptr_tensor + n_nodes*2);
    //auto end_time_copy = std::chrono::high_resolution_clock::now();
  if(model_type==torch::kFloat64) {
    positions =positions.transpose(0, 1).to(torch::kF64).contiguous();
  }
  else  if(model_type==torch::kFloat32) {
    positions =positions.transpose(0, 1).contiguous();
  }
    //auto end_time_trans = std::chrono::high_resolution_clock::now();
  //auto end_time_double = std::chrono::high_resolution_clock::now();
    #endif

  /*auto end_time_cont = std::chrono::high_resolution_clock::now();
  auto end_time_cout = std::chrono::high_resolution_clock::now();
  auto durationzz = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_zzytensor - start_time);
  std::cout << "创建指针与张量耗时: " << durationzz.count() << " milliseconds" << std::endl;
  auto durationcopy= std::chrono::duration_cast<std::chrono::milliseconds>(end_time_copy - end_time_zzytensor);
  std::cout << "转移环节耗时: " << durationcopy.count() << " milliseconds" << std::endl;
  auto durationtran = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_trans - end_time_copy);
  std::cout << "转置环节耗时: " << durationtran.count() << " milliseconds" << std::endl;
  auto durationdouble = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_double - end_time_trans);
  std::cout << "double化环节耗时: " << durationdouble.count() << " milliseconds" << std::endl;
  auto durationcont = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_cont - end_time_double);
  std::cout << "连续化环节耗时: " << durationcont.count() << " milliseconds" << std::endl;
  auto durationcout = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_cout - start_time);
  std::cout << "最后输出总耗时: " << durationcout.count() << " milliseconds" << std::endl;*/
}

void maceload::loadedges(){
    auto start_time = std::chrono::high_resolution_clock::now();
    // 1. 构建邻居列表
  //std::cout<<"bfbuild_positions[2][0]: " <<positions[2][0] <<"bfbuild_positions[2][0]device: "<<_mace_device_data->_d_px[2]<<std::endl;
    _neighbor_list = _neighbor_list_builder->Build();
  //std::cout<<"afbuild_positions[2][0]device: "<<_mace_device_data->_d_px[2]<<std::endl;
    auto end_build = std::chrono::high_resolution_clock::now();



    // 2. 计算总边数
  rbmd::Id n_edges = thrust::reduce(
    thrust::device,
    _neighbor_list->_d_neighbor_num.begin(),
    _neighbor_list->_d_neighbor_num.end(),
    rbmd::Id(0),
    thrust::plus<rbmd::Id>()
);
  std::cout <<"n_edges: "<< n_edges<<"\n";
//3


  thrust::device_vector<rbmd::Id> prefix_sum(n_nodes + 1);
/*  thrust::exclusive_scan(
      thrust::device,
      _neighbor_list->_d_neighbor_num.begin(),
      _neighbor_list->_d_neighbor_num.end(),
      prefix_sum.begin(),
      rbmd::Id(0)
  );
  thrust::device_vector<rbmd::Id> valid_indices(n_edges);
  auto valid_indices_ptr = thrust::raw_pointer_cast(valid_indices.data());
  auto start_idx_ptr = thrust::raw_pointer_cast(_neighbor_list->_start_idx.data());
  auto prefix_sum_ptr = thrust::raw_pointer_cast(prefix_sum.data());

  // 2. 使用 thrust::hip::par 策略
  thrust::for_each(
      thrust::hip::par,  // 关键：指定 HIP 后端
      thrust::make_counting_iterator(0),
      thrust::make_counting_iterator(n_edges),
      [=] __device__ (rbmd::Id idx) {
          auto it = thrust::upper_bound(
              thrust::hip::par,  // 统一执行策略
              prefix_sum_ptr,
              prefix_sum_ptr + n_nodes,
              idx
          );
          int node_id = it - prefix_sum_ptr - 1;
          valid_indices_ptr[idx] = start_idx_ptr[node_id] + (idx - prefix_sum_ptr[node_id]);
      }
  );*/


    // 3. 索引生成 - 使用循环
  thrust::device_vector<rbmd::Id> valid_indices(n_edges);

  // 获取原始指针
  //auto d_neighbor_num = _neighbor_list->_d_neighbor_num.data();
  auto d_start_idx = _neighbor_list->_start_idx.data();
  auto d_valid_indices = valid_indices.data();

  // 计算每个节点的写入偏移
  thrust::device_vector<rbmd::Id> write_offsets(n_nodes + 1);
  thrust::exclusive_scan(
      _neighbor_list->_d_neighbor_num.begin(),
      _neighbor_list->_d_neighbor_num.end(),
      write_offsets.begin()
  );
  write_offsets[n_nodes] = n_edges;
  auto d_write_offsets = write_offsets.data();

  // for_each,less edges
  thrust::for_each(
      thrust::make_counting_iterator(0),
      thrust::make_counting_iterator(n_nodes),
      [=] __device__ (rbmd::Id node) {
          rbmd::Id start = d_write_offsets[node];
          rbmd::Id end = d_write_offsets[node + 1];
          rbmd::Id start_val = d_start_idx[node];
          thrust::sequence(
              thrust::device,
              d_valid_indices + start,
              d_valid_indices + end,
              start_val
          );
      }
  );
  //transform,large edges
  /*thrust::transform(
    thrust::make_counting_iterator(0),
    thrust::make_counting_iterator(n_edges),
    valid_indices.begin(),
    [=] __device__ (rbmd::Id idx) {
        // 找到idx所属的节点
        rbmd::Id node = 0;
        while (idx >= d_write_offsets[node + 1]) node++;

        return d_start_idx[node] + (idx - d_write_offsets[node]);
    }
);*/
  //较慢的for循环
/*
  thrust::device_vector<rbmd::Id> valid_indices(n_edges);
  rbmd::Id startwrite = 0;
  rbmd::Id startwrite1,s_idx ;
  for (int i=0; i<n_nodes ;++i) {
    startwrite1 = startwrite+_neighbor_list->_d_neighbor_num[i];
    s_idx = _neighbor_list->_start_idx[i];
    thrust::sequence(valid_indices.begin()+startwrite, valid_indices.begin()+startwrite1, s_idx);
    startwrite = startwrite1;
  }*/



    auto end_indices = std::chrono::high_resolution_clock::now();

    // 5. 统一收集所有数据
    auto gather_data = [&](auto& dest, const auto& src) {
        dest.resize(n_edges);
        thrust::gather(
            thrust::device,
            valid_indices.begin(),
            valid_indices.end(),
            src.begin(),
            dest.begin()
        );
    };

    thrust::device_vector<rbmd::Id> edges1, edges2;
    thrust::device_vector<rbmd::Real> shiftx, shifty, shiftz;
    thrust::device_vector<rbmd::Real> unit_shiftx, unit_shifty, unit_shiftz;

    gather_data(edges1, _neighbor_list->_d_neighbors);
    gather_data(edges2, _neighbor_list->_d_neighbors_atoms);
    gather_data(shiftx, _neighbor_list->_d_shiftx);
    gather_data(shifty, _neighbor_list->_d_shifty);
    gather_data(shiftz, _neighbor_list->_d_shiftz);
    gather_data(unit_shiftx, _neighbor_list->_d_unit_shiftx);
    gather_data(unit_shifty, _neighbor_list->_d_unit_shifty);
    gather_data(unit_shiftz, _neighbor_list->_d_unit_shiftz);
    auto end_gather = std::chrono::high_resolution_clock::now();

    // 6. 优化张量创建 - 使用预处理指令
    auto model_type =r_max.dtype();
    #if USE_DOUBLE
        // 双精度路径

    auto temp_shifts = torch::empty({3, n_edges}, torch::kFloat64).to(torch::kCUDA);
    double* temp_shifts_data = temp_shifts.data_ptr<double>();
    thrust::copy(thrust::device, shiftx.begin(), shiftx.end(), temp_shifts_data);
    thrust::copy(thrust::device, shifty.begin(), shifty.end(), temp_shifts_data + n_edges);
    thrust::copy(thrust::device, shiftz.begin(), shiftz.end(), temp_shifts_data + 2 * n_edges);

    auto temp_unit_shifts = torch::empty({3, n_edges}, torch::kFloat64).to(torch::kCUDA);
    double* temp_unit_shifts_data = temp_unit_shifts.data_ptr<double>();
    thrust::copy(thrust::device, unit_shiftx.begin(), unit_shiftx.end(), temp_unit_shifts_data);
    thrust::copy(thrust::device, unit_shifty.begin(), unit_shifty.end(), temp_unit_shifts_data + n_edges);
    thrust::copy(thrust::device, unit_shiftz.begin(), unit_shiftz.end(), temp_unit_shifts_data + 2 * n_edges);
    if(model_type==torch::kFloat64) {
      shifts = temp_shifts.transpose(0, 1).contiguous();
      unit_shifts = temp_unit_shifts.transpose(0, 1).contiguous();
    }
    else if (model_type==torch::kFloat32) {
      shifts = temp_shifts.transpose(0, 1).to(torch::kFloat32).contiguous();
      unit_shifts = temp_unit_shifts.transpose(0, 1).to(torch::kFloat32).contiguous();
    }

        shifts = shifts.transpose(0, 1).contiguous();
        unit_shifts = unit_shifts.transpose(0, 1).contiguous();
    #else
        // 单精度路径
        // 创建临时单精度张量

     auto temp_shifts = torch::empty({3, n_edges}, torch::kFloat32).to(torch::kCUDA);
     float* temp_shifts_data = temp_shifts.data_ptr<float>();
     thrust::copy(thrust::device, shiftx.begin(), shiftx.end(), temp_shifts_data);
     thrust::copy(thrust::device, shifty.begin(), shifty.end(), temp_shifts_data + n_edges);
     thrust::copy(thrust::device, shiftz.begin(), shiftz.end(), temp_shifts_data + 2 * n_edges);

     auto temp_unit_shifts = torch::empty({3, n_edges}, torch::kFloat32).to(torch::kCUDA);
     float* temp_unit_shifts_data = temp_unit_shifts.data_ptr<float>();
     thrust::copy(thrust::device, unit_shiftx.begin(), unit_shiftx.end(), temp_unit_shifts_data);
     thrust::copy(thrust::device, unit_shifty.begin(), unit_shifty.end(), temp_unit_shifts_data + n_edges);
     thrust::copy(thrust::device, unit_shiftz.begin(), unit_shiftz.end(), temp_unit_shifts_data + 2 * n_edges);
  if(model_type==torch::kFloat64) {
     shifts = temp_shifts.transpose(0, 1).to(torch::kFloat64).contiguous();
     unit_shifts = temp_unit_shifts.transpose(0, 1).to(torch::kFloat64).contiguous();
   }
  else if (model_type==torch::kFloat32) {
    // 转置并转换为双精度
    shifts = temp_shifts.transpose(0, 1).contiguous();
    unit_shifts = temp_unit_shifts.transpose(0, 1).contiguous();

  }
    #endif

    // 7. 优化edge_index创建
    edge_index = torch::empty({2, n_edges}, torch::kInt64).to(torch::kCUDA);
    int64_t* edge_data = edge_index.data_ptr<int64_t>();

    #if USE_64BIT_IDS
        // 如果ID已经是64位，直接拷贝
        thrust::copy(thrust::device, edges1.begin(), edges1.end(), edge_data);
        thrust::copy(thrust::device, edges2.begin(), edges2.end(), edge_data + n_edges);
    #else
        // 如果ID是32位，需要转换
        thrust::transform(
            thrust::device,
            edges1.begin(),
            edges1.end(),
            edge_data,
            [] __device__ (rbmd::Id x) { return static_cast<int64_t>(x); }
        );

        thrust::transform(
            thrust::device,
            edges2.begin(),
            edges2.end(),
            edge_data + n_edges,
            [] __device__ (rbmd::Id x) { return static_cast<int64_t>(x); }
        );
       //edge_index = edge_index.to(torch::kInt32);
    #endif
    //auto end_tensors = std::chrono::high_resolution_clock::now();
  /*std::ofstream out_file("edges.txt");  // 创建并打开文件
if (out_file.is_open()) {
  out_file << edge_index << std::endl;  // 写入文件
  out_file.close();  // 关闭文件
} else {
  std::cerr << "无法打开文件！" << std::endl;
}*/

    // 8. 计时输出
    /*auto dur_build = std::chrono::duration_cast<std::chrono::milliseconds>(end_build - start_time);
    auto dur_indices = std::chrono::duration_cast<std::chrono::milliseconds>(end_indices - end_build);
    auto dur_gather = std::chrono::duration_cast<std::chrono::milliseconds>(end_gather - end_indices);
    auto dur_tensors = std::chrono::duration_cast<std::chrono::milliseconds>(end_tensors - end_gather);
    auto dur_total = std::chrono::duration_cast<std::chrono::milliseconds>(end_tensors - start_time);

    std::cout << "邻居构建耗时: " << dur_build.count() << " ms\n";
    std::cout << "计算边数，索引生成耗时: " << dur_indices.count() << " ms\n";
    std::cout << "数据收集耗时: " << dur_gather.count() << " ms\n";
    std::cout << "张量创建耗时: " << dur_tensors.count() << " ms\n";
    std::cout << "总耗时: " << dur_total.count() << " ms\n";*/

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
    //auto model_type =r_max.dtype();
    c10::Dict<std::string, torch::Tensor> input;
    //forces = torch::empty({ n_nodes, 3 }, torch_float_dtype.device(device));
    /*try
    {
      batch= batch.to(this->device);
      cell= cell.to(this->device);
      edge_index = edge_index.to(this->device);
      energy = energy.to(this->device);
      forces = forces.to(this->device);
      node_attrs = node_attrs.to(this->device);
      positions = positions.to(this->device);
      ptr=ptr.to(this->device);
      shifts = shifts.to(this->device);
      unit_shifts = unit_shifts.to(this->device);
      weight = weight.to(this->device);
    }
    catch (const c10::Error& e)
    {
      std::cerr << "Error in to device\n";
    }*/
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
    //auto start_time = std::chrono::high_resolution_clock::now();
    //model.eval();
    auto x_time = std::chrono::high_resolution_clock::now();
    //std::cout<<" positions[2][0]: "<< positions[2][0]<<" positions[2][0]device: "<<_mace_device_data->_d_px[2] <<" positions[2][0]device: "<<_mace_device_data->_d_py[2]<<" positions[2][0]device: "<<_mace_device_data->_d_pz[2]<<std::endl;

    //torch::NoGradGuard no_grad;
    c10::impl::GenericDict output = model.forward({ input, mask, vflag_global }).toGenericDict();
    auto y_time = std::chrono::high_resolution_clock::now();
    //auto end_time = std::chrono::high_resolution_clock::now();
    auto xytime = std::chrono::duration_cast<std::chrono::microseconds>(y_time - x_time);
    //auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
     std::cout << "xy: " << xytime.count() << " microseconds" << std::endl;
    //std::cout << "forward: " << duration.count() << " milliseconds" << std::endl;
    //std::cout << "模型前向完成" << std::endl;



    #if USE_DOUBLE
    energy = output.at("total_energy_local").toTensor();//gpu double
    double eng_vdwl = energy.item<double>();
    std::cout << eng_vdwl<< std::endl;
    forces = output.at("forces").toTensor();//gpu float
    auto forces_t = forces.transpose(0, 1).contiguous()*tokcalmol;
    double* forces_data_ptr = forces_t.data_ptr<double>();
    thrust::copy(forces_data_ptr, forces_data_ptr + n_nodes, _mace_device_data->_d_fx.begin());     // X方向
    thrust::copy(forces_data_ptr + n_nodes, forces_data_ptr + 2 * n_nodes, _mace_device_data->_d_fy.begin());  // Y方向
    thrust::copy(forces_data_ptr + 2 * n_nodes, forces_data_ptr + 3 * n_nodes, _mace_device_data->_d_fz.begin()); // Z方向
    if (vflag_global) 
    {
        vir = output.at("virials").toTensor();
        _mace_device_data->_d_virial_xx[0] += vir[0][0][0].item<double>();
        _mace_device_data->_d_virial_yy[0] += vir[0][1][1].item<double>();
        _mace_device_data->_d_virial_zz[0] += vir[0][2][2].item<double>();
        _mace_device_data->_d_virial_xy[0] += 0.5 * (vir[0][1][0].item<double>() + vir[0][0][1].item<double>());
        _mace_device_data->_d_virial_xz[0] += 0.5 * (vir[0][2][0].item<double>() + vir[0][0][2].item<double>());
        _mace_device_data->_d_virial_yz[0] += 0.5 * (vir[0][2][1].item<double>() + vir[0][1][2].item<double>());
    }
    #else
    energy = output.at("total_energy_local").toTensor().to(torch::kFloat);//gpu double
    //float eng_vdwl = energy.item<float>()*tokcalmol;
    float eng_vdwl = energy.item<float>();
    forces = output.at("forces").toTensor().to(torch::kFloat);//gpu float
    //auto forces_t = forces.transpose(0, 1).contiguous();//核心问题issue,爆炸
    auto forces_t = forces.transpose(0, 1).contiguous()*tokcalmol;
    float* forces_data_ptr = forces_t.data_ptr<float>();
    //float* forces_data_ptr = forces.data_ptr<float>();
    thrust::copy(forces_data_ptr, forces_data_ptr + n_nodes, _mace_device_data->_d_fx.begin());     // X方向
    thrust::copy(forces_data_ptr + n_nodes, forces_data_ptr + 2 * n_nodes, _mace_device_data->_d_fy.begin());  // Y方向
    thrust::copy(forces_data_ptr + 2 * n_nodes, forces_data_ptr + 3 * n_nodes, _mace_device_data->_d_fz.begin()); // Z方向
    //std::cout<<"force_x_yuan: "<<forces[2][0]<<"force_x"<<_mace_device_data->_d_fx[2]<<" positions[2][0]: "<< positions[2][0]<<" positions[2][0]device: "<<_mace_device_data->_d_px[2]<<std::endl;
    if (vflag_global)
    {
        vir = output.at("virials").toTensor().to(torch::kFloat);

        _mace_device_data->_d_virial[0] += vir[0][0][0].item<float>();
        _mace_device_data->_d_virial[1] += vir[0][1][1].item<float>();
        _mace_device_data->_d_virial[2] += vir[0][2][2].item<float>();
        _mace_device_data->_d_virial[3] += 0.5 * (vir[0][1][0].item<float>() + vir[0][0][1].item<float>());
        _mace_device_data->_d_virial[4] += 0.5 * (vir[0][2][0].item<float>() + vir[0][0][2].item<float>());
        _mace_device_data->_d_virial[5] += 0.5 * (vir[0][2][1].item<float>() + vir[0][1][2].item<float>());
    }
    #endif

  std::cout <<"energy:  "<<std::fixed << std::setprecision(6) <<eng_vdwl<< std::endl;
  auto interval = DataManager::getInstance().getConfigData()->Get<rbmd::Id>(
  "interval", "outputs", "thermo_out");

  std::ofstream outfile("thermo.txt", std::ios::app);
  if (outfile.tellp() == 0) {
    outfile << "step  mace_energy (kCal/mol)" << std::endl;
  }
  if (test_current_step % interval == 0) {
    outfile << test_current_step << "     " << eng_vdwl << std::endl;
  }
  outfile.close();
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
