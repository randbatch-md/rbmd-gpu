#pragma once
#include "object.h"

class Ensemble : public Object
{
public:
	Ensemble() {};
	virtual ~Ensemble()=default;

	virtual void Init() = 0;
	virtual int Run() = 0;

protected:
  
private:

};