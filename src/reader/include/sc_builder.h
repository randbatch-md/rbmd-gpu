#pragma once
#include "types.h"
#include "base_reader.h"
#include <memory>

class SCBuilder : public InBuilder
{
public:
	SCBuilder();
	virtual ~SCBuilder();

	int Build() override;
};