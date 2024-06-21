#include <iostream>
#include "command_line.h"
#include "version.h"

cxxopts::Options CommandLine::_opts("opts");

CommandLine::CommandLine(int argc, char* argv[])
{
	try
	{
		_co = _opts.parse(argc, argv);
	}
	catch (const std::exception&)
	{
		std::cout << "use \"-h\" to get correct command" << std::endl;
		return;
	}

	ParseCommand();
}

bool CommandLine::RunApplication()
{
	if (_co.count("j"))
	{
		auto file = _co["j"].as<std::string>();
		if (file != "rbmd.json")
		{
			_console->error("the json name must be rbmd.json!");
			return false;
		}

		return true;
	}

	return false;
}

void CommandLine::Initialize()
{
	_opts.add_options()("v", "version");
	_opts.add_options()("h", "help");
	_opts.add_options()("j", "json file", cxxopts::value<std::string>());
}

void CommandLine::ParseCommand()
{
	if (_co.count("v"))
	{
		std::cout << VERSION << std::endl;
	}
	else if (_co.count("h"))
	{
		std::cout << "-v: version" << std::endl;
		std::cout << "-j xxx.json: run rbmd" << std::endl;
	}
	else if (_co.count("j"))
	{
		return;
	}
	else
	{
		std::cout << "use \"-h\" to get more command" << std::endl;
	}
}
