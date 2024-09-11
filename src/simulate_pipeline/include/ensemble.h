#pragma once
#include "../common/object.h"
#include "velocity_controller.h"
#include "position_controller.h"
#include "force.h"
#include "shake_controller.h"
#include "temperature_controller.h"

class Ensemble : public Object
{
public:
	Ensemble() {};
	virtual ~Ensemble()=default;

	/**
	 *@brief init current ensemble
	*/
	virtual void Init() = 0;

	/**
	 *@brief Preprocessing settings before solving
	*/
	virtual void Presolve()=0;

	/**
	 * @brief Specific execution calculation part
	*/
	virtual void Solve() = 0;

	/**
	 * @brief Post processing after the completion of the current calculation step
	*/
	virtual void Postsolve() = 0;

	/**
	 * @brief The overall calculation process of the current calculation step
	 * @return Determine if there is an error here. If executed normally, return 0. If not, return other values
	*/
	int Run() 
	{
		Presolve();
		Solve();
		Postsolve();

		return 0;
	};

protected:

	std::shared_ptr<PositionController> _position_controller;
	std::shared_ptr<VelocityController> _velocity_controller;
	std::shared_ptr<Force> _force_controller;
	std::shared_ptr<ShakeController> _shake_controller;
	std::shared_ptr<TemperatureController> _temperature_controller;
};