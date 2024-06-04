#pragma once
#include "object.h"
#include "system.h"
#include "json/value.h"
#include "json/reader.h"
#include "Types.h"

class Executioner : public Object
{
public:
	Executioner(const Json::Value& node, std::shared_ptr<System>& system);
	virtual ~Executioner() = default;

protected:
	Json::Value _exec_node;
	std::shared_ptr<System>& _system;
	float _time_step;
	int _num_steps;
};