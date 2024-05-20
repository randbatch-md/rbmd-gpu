#pragma once
#include <string>

class BaseReader
{
public:
	BaseReader(const std::string& filePath) :
		_file_path(filePath){};

public:
	virtual int Execute() = 0;

protected:
	std::string _file_path;
};