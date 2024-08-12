#pragma once
#include "object.h"

class Ensemble : public Object
{
public:
	Ensemble() {};
	virtual ~Ensemble()=default;

	virtual void Init() = 0;
	virtual int Evolve() = 0;

protected:
	//std::unique_ptr<InitGlobal> _init_global;
  
private:

};