#pragma once
#include "vector"
#include "types.h"
class ForceController
{
public:
	ForceController() {};
	virtual ~ForceController()=default;

	virtual void Update()=0;
	virtual void Init() {};
	virtual int Execute() = 0;

public:
	std::vector<rbmd::Real> _fx;
	std::vector<rbmd::Real> _fy;
	std::vector<rbmd::Real> _fz;
};