#pragma once

#include <memory>
#include "object.h"

class CommandLine;
class BaseReader;
class System;

class Application : public Object
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
	std::shared_ptr<System> _system;

};