// #pragma once
// #include "structure_reader.h"
//
// class Charge_Reader : public StructureReder
//{
// public:
//	Charge_Reader(const std::string& filePath, MDData& data);
//	~Charge_Reader() = default;
//
// protected:
//	int ReadData() override;
//	void AllocateDataSpace() override;
//
// private:
//	int ReadAtoms(const rbmd::Id& atoms_num);
//	int ReadVelocity(const rbmd::Id& atoms_num);
// };