#pragma once
#include "ensemble.h"

class NVPensemble : public Ensemble
{
public:
	NVPensemble();
	virtual~NVPensemble() = default;

protected:
	void Init() override;
	void Presolve() override;
	void Solve() override;
	void Postsolve() override;

private:

};