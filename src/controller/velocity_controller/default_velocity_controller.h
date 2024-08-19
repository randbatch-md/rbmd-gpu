#pragma once
#include "velocity_controller.h"

class DefaultVelocityController : VelocityController
{
public:
	DefaultVelocityController();
	virtual ~DefaultVelocityController()=default;

	void Update() override;
	void Init() override;
private:
	rbmd::Id _num_atoms;
	rbmd::Real _dt;
	rbmd::Real _fmt2v;//这个参数流的一个接口 原本参数是一个结构体 _unit_factor._fmt2v	
};