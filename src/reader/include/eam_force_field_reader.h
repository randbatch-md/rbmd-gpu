#pragma once

#include "../include/force_field_reader.h"

class EAMForceFieldReader : public ForceFieldReader
{
public:
	EAMForceFieldReader(MDData& data);
	~EAMForceFieldReader() = default;

protected:
	int ReadData() override;
	void AllocateDataSpace() override;

private:
	int ReadMass(const rbmd::Id& numAtomTypes);
	int ReadEAM();

};