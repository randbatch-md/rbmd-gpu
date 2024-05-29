#pragma once

#include <memory>

class CommandLine;
class BaseReader;
class System;

class Application
{
public:
	Application(int argc, char* argv[]);
	virtual ~Application() = default;

public:
	int Run();

protected:
	virtual int Execute() = 0;

protected:
	std::shared_ptr<CommandLine> _command_line;
	std::shared_ptr<BaseReader> _reader;
	std::shared_ptr<System> _system;

};