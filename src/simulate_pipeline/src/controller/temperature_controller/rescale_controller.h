#pragma once
#include "temperature_controller.h"

class RescaleController : public TemperatureController
{
public:
	RescaleController();
	virtual ~RescaleController()=default;

	void Init() override;
	void Update() override;
	void ComputeTemp();
	void UpdataVelocity();
private:
	rbmd::Id _num_atoms;
	rbmd::Real _mvv2e;
	rbmd::Real _kB;
	rbmd::Real _temp_sum;
	rbmd::Real _temp;
};