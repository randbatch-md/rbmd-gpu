#pragma once
#include "ensemble.h"

class NVPensemble : public Ensemble
{
public:
	NVPensemble();
	virtual~NVPensemble() = default;

protected:
	void Init() override;
	int Run() override;

private:

};