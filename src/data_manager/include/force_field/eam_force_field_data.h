#pragma once

#include <thrust/host_vector.h>

#include <fstream>

#include "common/types.h"
#include "force_field_data.h"

class EAMForceFieldData : public ForceFieldData {
 public:
  bool checkForceField() const override { return true; }

  void ReadPotentialFile(const std::string& filename)
  {
    std::ifstream input_file(filename);
    if (!input_file.is_open())
    {
      std::cerr << "Unable to open the file." << std::endl;
    }
    for (int i = 0; i < 2; ++i)
    {
      input_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    input_file >> _nrho >> _drho >> _nr>> _dr >> _cut_off;

    _h_frho.resize(_nrho + 1);
    _h_zr.resize(_nr+ 1);
   _h_rhor.resize(_nrho + 1);

    for (int i = 0; i < _nrho; ++i)
    {
      input_file >> _h_frho[i];
    }

    for (int i = 0; i < _nr; ++i)
    {
      input_file >> _h_zr[i];
    }

    for (int i = 0; i < _nrho; ++i)
    {
      input_file >> _h_rhor[i];
    }
    input_file.close();
  }

  rbmd::Real _cut_off;

  /// F(ρ) on host
  rbmd::Id _nrho;
  rbmd::Real _drho;
  //rbmd::Real* _h_frho;
  thrust::host_vector<rbmd::Real> _h_frho;


  /// ρ(r) on host
  // rbmd::Real* _h_rhor;
  thrust::host_vector<rbmd::Real> _h_rhor;

  /// ϕ(r) on host
  rbmd::Id _nr;
  rbmd::Real _dr;
  // rbmd::Real* _h_zr;
  thrust::host_vector<rbmd::Real> _h_zr;

};
