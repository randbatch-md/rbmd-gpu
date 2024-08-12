#pragma once
#include "ensemble.h"

class NVTensemble : public Ensemble
{
public:
	NVTensemble();
	~NVTensemble();

protected:
	void Init();
	int Execute();

private:

};