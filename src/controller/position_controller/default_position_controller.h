#pragma once
#include "position_controller.h"

class DefaultPositionController : PositionController
{
public:
	DefaultPositionController();
	virtual ~DefaultPositionController()=default;

protected:
	int Update() override;
private:

};