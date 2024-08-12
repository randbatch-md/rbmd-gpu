#pragma once
#include "ensemble.h"

class NVPensemble : public Ensemble
{
public:
	NVPensemble();
	~NVPensemble();

protected:
	void Init();
	int Execute();

private:

};