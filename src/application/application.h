#pragma once

#include <memory>

class CommandLine;
class BaseReader;
class System;

class Application
{
public:
	virtual void Run() = 0;

protected:
	std::shared_ptr<CommandLine> _command_line;
	std::shared_ptr<BaseReader> _reader;
	std::shared_ptr<System> _system;

};