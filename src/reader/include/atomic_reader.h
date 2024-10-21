#pragma once
#include "../common/types.h"
#include "structure_reader.h"
#include <map>

class AtomicReader : public StructureReder {
 public:
  AtomicReader(const std::string& filePath, MDData& data);
  ~AtomicReader() = default;

 protected:
  int ReadData() override;
  void AllocateDataSpace() override;

 private:
  int ReadAtoms(const rbmd::Id& atoms_num);
  int ReadVelocity(const rbmd::Id& atoms_num);
  int ReadBond(const rbmd::Id& atoms_num);
  int ReadAngle(const rbmd::Id& atoms_num);
  int ReadDihedrals(const rbmd::Id& atoms_num);
  void SetSpecialBonds();

  std::multimap<rbmd::Id, rbmd::Id> _special_map;
};