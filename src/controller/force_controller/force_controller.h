#pragma once

class ForceController
{
public:
	ForceController() {};
	virtual ~ForceController()=default;

protected:
	virtual int Init() {};

	virtual int Update() {};
	virtual int Update()=0;


};