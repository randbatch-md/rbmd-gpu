#pragma once
#include "object.h"
#include "velocity_controller.h"
#include "position_controller.h"
#include "force_controller.h"

class Ensemble : public Object
{
public:
	Ensemble() {};
	virtual ~Ensemble()=default;

	virtual void Init() = 0;
	virtual int Run() = 0;

protected:

	std::shared_ptr<PositionController> _position_controller;
	std::shared_ptr<VelocityController> _velocity_controller;
	std::shared_ptr<ForceController> _force_controller;


};