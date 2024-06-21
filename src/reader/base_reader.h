#pragma once
#include <string>
#include "object.h"

class BaseReader : public Object
{
public:
	BaseReader(const std::string& filePath) :
		_file_path(filePath){};

public:
	virtual int Execute() = 0;
	void SetFile(const std::string& file) { _file_path = file; }

protected:
	std::string _file_path;
};