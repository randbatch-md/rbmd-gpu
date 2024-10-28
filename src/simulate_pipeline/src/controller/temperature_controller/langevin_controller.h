#pragma once
#include "temperature_controller.h"
#include "simulate.h"

class LangevinController : public TemperatureController
{
public:
	LangevinController();
	virtual ~LangevinController();

	void Init() override;
	void Update() override;

	/**
	 * @brief Calculate the current stage temperature
	*/
	void ComputeTemp() override;

	/**
	 * @brief Update current speed through temperature
	*/
	void UpdataForce();

	rbmd::Real FetchSample_1D(const bool& random, const rbmd::Real& mu, const rbmd::Real& sigma);

	rbmd::Real RandomValue(const rbmd::Real& Min, const rbmd::Real& Max);

private:
	rbmd::Real _mvv2e;
	rbmd::Real _kB;
    rbmd::Real* _d_temp_contrib;
	rbmd::Real _dt;
};