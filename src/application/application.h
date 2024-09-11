#pragma once

#include <memory>
#include "../common/object.h"
#include "../data_manager/include/config_data.h"
#include "ensemble.h"
class CommandLine;
class BaseReader;
//class System;
//class JsonParser;
class Executioner;
class Ensemble;

class Application : public Object
{
public:
	/**
	 * @brief 
	 * @param argc 
	 * @param argv 
	*/
	Application(int argc, char* argv[]);
	virtual ~Application() = default;

public:
	int Run();

protected:
	virtual int Execute() = 0;

private:
	bool Check();

protected:
	//std::shared_ptr<CommandLine> _command_line;
	//std::shared_ptr<System> _system;
    std::vector<std::shared_ptr<Ensemble>> _simulate_pipelines; 
	std::shared_ptr<Ensemble> _simulate_pipeline;

	std::vector<Json::Value> _simulate_nodes;

	std::shared_ptr<ConfigData> _config_data;

	std::shared_ptr<Executioner> _executioner;
};