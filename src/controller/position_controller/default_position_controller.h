#pragma once
#include "position_controller.h"

class DefaultPositionController : PositionController
{
public:
	DefaultPositionController();
	virtual ~DefaultPositionController()=default;

	void Update() override;
	void Init() override;

private:

};