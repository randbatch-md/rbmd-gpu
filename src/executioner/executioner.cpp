#include "executioner.h"

Executioner::Executioner(const Json::Value& node, std::shared_ptr<System>& system) :
	_exec_node(node),
	_system(system),
	_time_step(_exec_node["timestep"].asFloat()),
	_num_steps(_exec_node["num_steps"].asInt())
{

}

void Executioner::Init()
{

}

int Executioner::Execute()
{
	_system->Evolve();

	return 0;
}
