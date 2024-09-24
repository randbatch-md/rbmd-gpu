#include "neighbor_list.h"

#include <thrust/host_vector.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "common/rbmd_define.h"
#include "data_manager.h"
#include "linked_cell_locator.h"

NeighborList::NeighborList(rbmd::Id total_atoms_num, bool is_half) {
  CHECK_RUNTIME(
      MALLOCHOST(reinterpret_cast<void **>(&(this->_h_is_half)), sizeof(bool)));
  *(this->_h_is_half) = is_half;
  this->_d_neighbor_num.resize(total_atoms_num);
  this->_d_max_neighbor_num.resize(total_atoms_num);
  this->_start_idx.resize(total_atoms_num);
  this->_end_idx.resize(total_atoms_num);
  CHECK_RUNTIME(MALLOC(&this->_d_is_half, sizeof(bool)));
  CHECK_RUNTIME(MEMCPY(this->_d_is_half, this->_h_is_half, sizeof(bool), H2D));
}

void NeighborList::print(const std::string& filename) {
  // 准备需要的数据，假设这些是你已经准备好的数据
  thrust::host_vector<rbmd::Id> atomid2idx = LinkedCellLocator::GetInstance().GetLinkedCell()->_atom_id_to_idx;
  thrust::host_vector<rbmd::Id> atoms_id = DataManager::getInstance().getDeviceData()->_d_atoms_id;
  thrust::host_vector<rbmd::Id> start_idx = _start_idx;
  thrust::host_vector<rbmd::Id> end_idx = _end_idx;
  thrust::host_vector<rbmd::Id> _neighbors = _d_neighbors;

  // 打开文件流
  std::ofstream outFile(filename);
  if (!outFile.is_open()) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return;
  }

  // 计算最大邻居数量
  int max_neighbors = 0;
  for (int i = 0; i < atomid2idx.size(); ++i) {
    int num_neighbors = end_idx[atomid2idx[i]] - start_idx[atomid2idx[i]];
    if (num_neighbors > max_neighbors) {
      max_neighbors = num_neighbors;
    }
  }

  // 输出 CSV 表头
  for (int i = 0; i < atomid2idx.size(); ++i) {
    outFile << "neighbor_" << i;
    if (i != atomid2idx.size() - 1) {
      outFile << ",";
    }
  }
  outFile << std::endl;

  // 输出 CSV 每一行，处理不同邻居数量
  for (int row = 0; row < max_neighbors; ++row) {
    for (int i = 0; i < atomid2idx.size(); ++i) {
      int start = start_idx[atomid2idx[i]];
      int end = end_idx[atomid2idx[i]];
      int num_neighbors = end - start;

      // 如果当前原子有足够的邻居，输出邻居ID，否则输出空值
      if (row < num_neighbors) {
        outFile << atoms_id[_neighbors[start + row]];
      } else {
        outFile << "";  // 空值
      }

      if (i != atomid2idx.size() - 1) {
        outFile << ","; // 逗号分隔符
      }
    }
    outFile << std::endl; // 行结束
  }

  // 关闭文件流
  outFile.close();
}
