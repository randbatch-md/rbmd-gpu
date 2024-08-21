//#pragma once
//#include "structure_reader.h"
//
//class FullReader : public StructureReder
//{
//public:
//	FullReader(const std::string& filePath, MDData& data);
//	~FullReader() = default;
//
//protected:
//	int ReadData() override;
//	void AllocateDataSpace() override;
//
//private:
//	int ReadAtoms(const rbmd::Id& atoms_num);
//	int ReadVelocity(const rbmd::Id& atoms_num);
//};