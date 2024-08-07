#pragma once
#include <memory>

struct StructureData;
struct StructureInfoData;
struct ForceFieldData;
struct LJStructureData;

struct MDData
{
	MDData()
	{
		_structure_data = std::make_shared<LJStructureData>();
		_structure_info_data = std::make_shared<StructureInfoData>();
		_force_field_data = std::make_shared<ForceFieldData>();
	}

	std::shared_ptr<StructureData> _structure_data;
	std::shared_ptr<StructureInfoData> _structure_info_data;
	std::shared_ptr<ForceFieldData> _force_field_data;
};
