#include "application.h"
#include "command_line.h"

Application::Application(int argc, char* argv[])
{
	CommandLine::Initialize();
	_command_line = std::make_shared<CommandLine>(argc, argv);
}

int Application::Run()
{
	if (_command_line->RunApplication())
	{
		return Execute();
	}

	return 0;
}
