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

  void MolecularMapInsert(rbmd::Id& key, rbmd::Id& value);
  void AtomstoMolecular(rbmd::Id& key, rbmd::Id& value);
  void SetMolecularGroup();

  std::multimap<rbmd::Id, rbmd::Id> _special_map;
  std::map<rbmd::Id, std::vector<rbmd::Id>> _molecular_map;
  std::map<rbmd::Id, rbmd::Id> atom_to_molecular_map;
};