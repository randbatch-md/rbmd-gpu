#pragma once
#include "ensemble.h"

class NVTensemble : public Ensemble
{
public:
	NVTensemble();
	virtual ~NVTensemble() = default;

protected:
	void Init() override;
	int Run() override;

private:

};