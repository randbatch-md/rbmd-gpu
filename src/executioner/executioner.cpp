#include "executioner.h"

Executioner::Executioner(const Json::Value& node, std::vector<std::shared_ptr<Ensemble>>& simulate_pipeline) :
	_exec_node(node),
	//_system(system),
	_simulate_pipeline(simulate_pipeline)
	_time_step(_exec_node["timestep"].asFloat()),
	_num_steps(_exec_node["num_steps"].asInt())
{

}

void Executioner::Init()
{

}

int Executioner::Execute()
{
	_system->Evolve(); //dynamic

	return 0;
}
