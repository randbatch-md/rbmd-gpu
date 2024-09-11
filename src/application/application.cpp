#include "application.h"
#include "command_line.h"
//#include "JsonParser.h"
#include "data_manager/include/data_manager.h"

Application::Application(int argc, char* argv[])
{
	CommandLine::Initialize();
	//_command_line = std::make_shared<CommandLine>(argc, argv);
	_config_data = DataManager::getInstance().getConfigData(); 
}

int Application::Run()
{
	//if (_command_line->RunApplication())
	//{
	//	//_config_data = std::make_shared<ConfigData>("rbmd.json");
	//	if (!Check())
	//	{
	//		_console->error("the data name must be rbmd.data!");
	//		return -1;
	//	}
	//
	//	return Execute();
	//}

	return Execute();
	return 0;
}

bool Application::Check()
{
	auto file = _config_data->Get<std::string>("file", "init_configuration", "read_data");

	if (file != "rbmd.data")
	{
		return false;
	}

	return true;
}
