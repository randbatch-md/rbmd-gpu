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
	rbmd::Real _fmt2v;//�����������һ���ӿ� ԭ��������һ���ṹ�� _unit_factor._fmt2v	
};