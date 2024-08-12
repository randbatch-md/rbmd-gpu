#pragma once
#include "object.h"

class Ensemble : public Object
{
public:
	Ensemble();
	~Ensemble();

	virtual int Evolve() = 0;

protected:
	std::unique_ptr<InitGlobal> _init_global;
  
private:

};