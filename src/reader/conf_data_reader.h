#pragma once
#include "mmap_reader.h"
#include "model/data_header.h"
#include "model/data_body.h"

class ConfDataReder : public MmapReader
{
public:
	ConfDataReder(const std::string& filePath);
	~ConfDataReder() = default;

	int Execute() override;

public:
	auto& GetData() { return _Data; }

protected:
	virtual int ReadData() = 0;

private:
	int ReadHeader();

private:
	std::unique_ptr<DataHeader> _header;
	std::unique_ptr<DataBody> _Data;
};
