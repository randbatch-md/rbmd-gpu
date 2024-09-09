#include "executioner.h"

Executioner::Executioner(const Json::Value& node, std::shared_ptr<Ensemble>& simulate_pipeline) :
	_exec_node(node),
	//_system(system),
	_simulate_pipeline(simulate_pipeline),
	_time_step(_exec_node["num_steps"].asFloat()),
	_num_steps(_exec_node["num_steps"].asInt())
{
	_current_step = 0;
	_current_time = 0;
}

void Executioner::Init()
{
	_simulate_pipeline->Init();
}

int Executioner::Execute()
{
	while (KeepGoing())
	{
		_current_step++;
		_current_time += _time_step;
		//_system->Evolve(); //dynamic
		_simulate_pipeline->Run();
		return 0;
	}
}

bool Executioner::KeepGoing()
{
	bool keep_going = false;

	if (_current_step <= _num_steps)
		keep_going = true;

	return keep_going;
}
