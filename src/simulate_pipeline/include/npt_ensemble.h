#pragma once
#include "ensemble.h"

class NPTensemble : public Ensemble
{
public:
	NPTensemble();
	virtual~NPTensemble() = default;

protected:
	void Init() override;
	void Presolve() override;
	void Solve() override;
	void Postsolve() override;

private:

};