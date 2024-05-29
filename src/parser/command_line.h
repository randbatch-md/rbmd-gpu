#pragma once

#include "cxxopts.hpp"

class CommandLine
{
public:
	CommandLine(int argc, char* argv[]);
	virtual ~CommandLine() = default;

public:
	bool RunApplication();
	static void Initialize();
	std::string GetFile();

private:
	void ParseCommand();

private:
	cxxopts::ParseResult _co;
	static cxxopts::Options _opts;

};