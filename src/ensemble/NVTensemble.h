#pragma once
#include "ensemble.h"

class NVTensemble : public Ensemble
{
public:
	NVTensemble();
	virtual ~NVTensemble() = default;

protected:
	void Init() override;
	void Presolve() override;
	void Solve() override;
	void Postsolve() override;

private:


};