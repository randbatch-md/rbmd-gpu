#pragma once
#include "position_controller.h"

class DefaultPositionController : PositionController
{
public:
	DefaultPositionController();
	virtual ~DefaultPositionController()=default;

	void Init() override;
	void Update() override;

	void SetCenterTargetPositions();

private:
	rbmd::Id _num_atoms;
	rbmd::Real _dt;
	std::string _init_type;
};