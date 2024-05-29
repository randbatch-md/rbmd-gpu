#pragma once
#include "application.h"

class MDApplication : public Application
{
public:
	MDApplication(int argc, char* argv[]);
	~MDApplication() = default;

public:
	int Execute () override;

};