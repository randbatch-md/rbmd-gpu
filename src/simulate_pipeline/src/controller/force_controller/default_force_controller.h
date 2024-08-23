#pragma once
#include "force_controller.h"

class DefaultForceController : public ForceController
{
public:
	DefaultForceController();
	virtual ~DefaultForceController()=default;

protected:
	void Update() override;
	void Init() override;
	int Execute() override;


private:

};