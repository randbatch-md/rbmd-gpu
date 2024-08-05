#pragma once
#include "object.h"
#include "system.h"
#include "json/value.h"
#include "json/reader.h"
#include "types.h"

class Executioner : public Object
{
public:
	Executioner(const Json::Value& node, std::shared_ptr<System>& system);
	virtual ~Executioner() = default;

public:
	void Init();
	int Execute();

protected:
	Json::Value _exec_node;
	std::shared_ptr<System>& _system;
	float _time_step;
	int _num_steps;
};