#pragma once

#include <fcntl.h>
#include <string>

class BaseReader
{
public:
	BaseReader(const std::string& filePath) :
	file_path(filePath)
	{

	}

public:
	virtual void Execute() = 0;

protected:
	std::string file_path;
};