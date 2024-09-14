#pragma once
#include "velocity_controller.h"
#include "../executioner/executioner.h"
class DefaultVelocityController : public VelocityController
{
public:
	DefaultVelocityController();
	virtual ~DefaultVelocityController()=default;

	void Init() override;
	void Update() override;
	void Update2() override;

private:
	rbmd::Id _num_atoms;
	rbmd::Real _dt;
	rbmd::Real _fmt2v;
};