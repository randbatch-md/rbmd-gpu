#pragma once
#include "ensemble.h"

class NVEensemble : public Ensemble
{
public:
	NVEensemble();
	virtual ~NVEensemble() = default;

protected:
	void Init() override;
	int Evolve() override;
private:

};