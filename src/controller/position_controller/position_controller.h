#pragma once

class PositionController
{
public:
	PositionController() {};
	virtual ~PositionController()=default;

protected:
	virtual int Update() {};

};