#pragma once
 
#include <c10/util/Logging.h>
#include <algorithm>
#include <iostream>
#include <array>
#include <vector>
#include "force.h"
#include <chrono>
#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include "../../data_manager/include/data_manager.h"
#include "../../data_manager/include/model/device_data.h"
#include "../../data_manager/include/model/box.h"
#include "../../data_manager/include/scheduler/memory_scheduler.h"
#include "../../data_manager/include/model/structure_info_data.h"
#include "../neighbor_list/include/linked_cell/linked_cell.h"
#include "neighbor_list/include/neighbor_list/neighbor_list.h"
#define GLOG_USE_GLOG_EXPORT
#pragma push_macro("INFO")
#pragma push_macro("WARNING")
#pragma push_macro("ERROR")
#include <glog/logging.h>
#pragma pop_macro("INFO")
#pragma pop_macro("WARNING")
#pragma pop_macro("ERROR")
#include "neighbor_list/include/neighbor_list_builder/mace_neighbor_list_builder.h"
#include <torch/torch.h>
#undef REDUCE
#include <torch/script.h>
#define REDUCE hipcub::DeviceReduce::Sum


class maceload  : public Force
{
    public:
        maceload();
        virtual ~maceload()=default;
        //maceload(){};
        void Init() override;
        void loadpositions(); 
        void loadcell();
        //void loadatoms();
        void loadatoms(std::vector<std::string>);
        void loadmass(thrust::device_vector<rbmd::Real>);
        void loadedges();
        void loadedges1();
        void load_atomid();
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
            rbmd::Id* neighbor_cell, rbmd::Id neighbor_cell_num);
        void  Execute() override;
        void forward();
        double energyout(c10::impl::GenericDict ouput);
        torch::Tensor virialout(c10::impl::GenericDict ouput);
        torch::Tensor atomenergy(c10::impl::GenericDict ouput);
        torch::Tensor vir;
        c10::TensorOptions torch_float_dtype = torch::dtype(torch::kF64);
        c10::TensorOptions torch_int_dtype = torch::dtype(torch::kInt64);
        double mace_r_max;
        torch::jit::script::Module model;
        int n_nodes, n_edges;
        bool vflag_global=false;
        float tokcalmol = 23.0605;
        //double eng_vdwl;
        //ContPointLocator _locator;
    private:
        std::shared_ptr<LinkedCell> _linked_cell;
        c10::Device device = torch::Device(torch::DeviceType::CUDA);
        //caffe2::TypeMeta model_type =double;
        torch::Tensor r_max;
        std::vector<int> atoms_table;
        torch::Tensor d_atoms_table;
        std::vector<std::string> mace_feats_table;
        std::vector<float> mace_feats_mass_table;//thrust::device_vector<float>这个类型
        int64_t n_node_feats;
        torch::Tensor cell;
        torch::Tensor weight;
        torch::TensorOptions options;
        torch::Tensor energy;
        torch::Tensor positions;
        torch::Tensor positionx,positiony,positionz;
        torch::TensorOptions optsint32 = torch::TensorOptions().dtype(torch::kInt32);
        torch::TensorOptions optsint64 = torch::TensorOptions().dtype(torch::kInt64);
        torch::TensorOptions optsfloat32 = torch::TensorOptions().dtype(torch::kF32);
        torch::TensorOptions optsfloat64 = torch::TensorOptions().dtype(torch::kF64);
        //using DeviceType = Kokkos::Experimental::HIP;
        torch::Tensor batch;
        torch::Tensor forces;
        torch::Tensor ptr = torch::zeros({ 2 }, torch::dtype(torch::kInt64));
        torch::Tensor node_attrs;
        torch::Tensor edge_index;
        torch::Tensor unit_shifts;
        torch::Tensor shifts;
        torch::Tensor mask;
        double eng_vdwl;
        const std::array<std::string, 118> periodic_table =
        { "H", "He",
         "Li", "Be",                                                              "B",  "C",  "N",  "O",  "F", "Ne",
         "Na", "Mg",                                                             "Al", "Si",  "P",  "S", "Cl", "Ar",
         "K",  "Ca", "Sc", "Ti",  "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
         "Rb", "Sr",  "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te",  "I", "Xe",
         "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
                           "Hf", "Ta",  "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
         "Fr", "Ra", "Ac", "Th", "Pa",  "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
                           "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og" };
        const std::array<float, 118> mass_periodic_table =
        { 1.008, 4.0026,
         6.94  , 9.0122,                                                                                  10.81,  12.011, 14.007,  15.999,  18.998, 20.180,
         22.990, 24.305,                                                                                  26.982, 28.085, 30.974,  32.06,   35.45,  39.948,
         39.098, 40.078, 44.956, 47.867,  50.942, 51.996,  54.938, 55.845, 58.933, 58.693, 63.546, 65.38, 69.723, 72.630, 74.922,  78.971,  79.904, 83.798,
         85.468, 87.62,  88.906, 91.224,  92.906, 95.95,   98,     101.07, 102.91, 106.42, 107.87, 112.41, 114.82, 118.71, 121.76, 127.60,  126.90, 131.29,
         132.91, 137.33, 138.91, 140.12,  140.91, 144.24,  145,    150.36, 151.96, 157.25, 158.93, 162.50, 164.93, 167.26, 168.93, 173.05,  174.97,
                                 178.49, 180.95,  183.84,  186.21, 190.23, 192.22, 195.08, 196.97, 200.59, 204.38, 207.2,  208.98,  209,    210,     222,
         223,    226,    227,    232.04,  231.04, 238.03,  237,    244,    243,    247,    247,    251,    252,    257,    258,     259,    266,
                                 267,    268,     269,     270,    277,    278,    281,    282,    285,    286,    289,    290,     293,    294,     294 };
        /*this->_structure_info_data =
        DataManager::getInstance().getMDData()->_structure_info_data;
        this->_device_data = DataManager::getInstance().getDeviceData();*/
        std::shared_ptr<DeviceData> _mace_device_data;
        std::shared_ptr<NeighborList> _neighbor_list = nullptr;
        //std::shared_ptr<StructureData> _structure_data;
        //std::shared_ptr<StructureInfoData> _mace_structure_info_data;
        std::shared_ptr<BaseNeighborListBuilder> _neighbor_list_builder;
        std::shared_ptr<LinkedCell> _mace_linked_cell;
        std::shared_ptr<NeighborList> _list;

  struct GenerateRange {
    thrust::device_vector<rbmd::Id>::const_iterator a_begin;
    thrust::device_vector<rbmd::Id>::const_iterator b_begin;

    GenerateRange(thrust::device_vector<rbmd::Id>::const_iterator a_begin, thrust::device_vector<rbmd::Id>::const_iterator b_begin)
        : a_begin(a_begin), b_begin(b_begin) {}

    __host__ __device__
    int operator()(const rbmd::Id& idx) const {
      return *(a_begin + idx);
    }
  };

  struct EndRange {
    thrust::device_vector<rbmd::Id>::const_iterator b_begin;

    EndRange(thrust::device_vector<rbmd::Id>::const_iterator b_begin) : b_begin(b_begin) {}

    __host__ __device__
    bool operator()(const rbmd::Id& val, const rbmd::Id& idx) const {
      return val >= *(b_begin + idx);
    }
  };

};

extern maceload macetest;
