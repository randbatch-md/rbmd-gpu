#include "application.h"
#include "command_line.h"
#include "JsonParser.h"

Application::Application(int argc, char* argv[])
{
	CommandLine::Initialize();
	_command_line = std::make_shared<CommandLine>(argc, argv);
}

int Application::Run()
{
	if (_command_line->RunApplication())
	{
		_parser = std::make_shared<JsonParser>("rbmd.json");
		if (!Check())
		{
			_console->error("the data name must be rbmd.data!");
			return -1;
		}

		return Execute();
	}

	return 0;
}

bool Application::Check()
{
	auto file = _parser->Get<std::string>("file", "init_configuration", "read_data");

	if (file != "rbmd.data")
	{
		return false;
	}

	return true;
}
