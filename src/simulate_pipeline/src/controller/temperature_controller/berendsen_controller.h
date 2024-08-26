#pragma once
#include "temperature_controller.h"

class BerendsenController : public TemperatureController
{
public:
	BerendsenController();
	virtual ~BerendsenController()=default;

	void Init() override;
	void Update() override;
	void ComputeTemp();
	void UpdataVelocity();
private:
	rbmd::Id _num_atoms;
	rbmd::Real _dt;
	rbmd::Real _mvv2e;
	rbmd::Real _kB;
	rbmd::Real _temp_sum;
	rbmd::Real _temp;

	rbmd::Real _Tdamp;
};