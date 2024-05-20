#pragma once
#include "conf_data_reader.h"

class AtomicReader : public ConfDataReder
{
public:
	AtomicReader(const std::string& filePath);
	~AtomicReader() = default;

protected:
	int ReadData() override;
};