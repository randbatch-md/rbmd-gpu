#pragma once

#include "cxxopts.hpp"
#include "object.h"

class CommandLine : public Object
{
public:
	CommandLine(int argc, char* argv[]);
	virtual ~CommandLine() = default;

public:
	bool RunApplication();
	static void Initialize();

private:
	void ParseCommand();

private:
	cxxopts::ParseResult _co;
	static cxxopts::Options _opts;

};