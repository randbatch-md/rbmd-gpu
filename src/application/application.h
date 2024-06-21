#pragma once

#include <memory>
#include "object.h"

class CommandLine;
class BaseReader;
class System;
class JsonParser;
class Executioner;

class Application : public Object
{
public:
	Application(int argc, char* argv[]);
	virtual ~Application() = default;

public:
	int Run();

protected:
	virtual int Execute() = 0;

private:
	bool Check();

protected:
	std::shared_ptr<CommandLine> _command_line;
	std::shared_ptr<System> _system;
	std::shared_ptr<JsonParser> _parser;
	std::shared_ptr<Executioner> _executioner;
};