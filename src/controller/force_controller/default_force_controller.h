#pragma once
#include "force_controller.h"

class DefaultForceController : ForceController
{
public:
	DefaultForceController();
	virtual ~DefaultForceController()=default;

protected:
	int Update() override;
	int Init() override;
	int Execute() override;


private:

};