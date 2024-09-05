#pragma once
#include "types.h"
#include "data_builder.h"
#include <memory>

class FileBuilder : public DataBuilder
{
public:
	FileBuilder(const std::string& filePath, const std::string& atom_style, const std::string& force_field);
	virtual ~FileBuilder();

	int Build() override;

protected:
	std::shared_ptr<BaseReader> _reader;
};