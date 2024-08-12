#pragma once
#include "velocity_controller.h"

class DefaultVelocityController : VelocityController
{
public:
	DefaultVelocityController();
	virtual ~DefaultVelocityController()=default;

protected:
	int Update() override;
private:

};