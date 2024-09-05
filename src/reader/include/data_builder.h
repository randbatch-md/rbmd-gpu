#pragma once
#include <string>
#include "object.h"

class DataBuilder : public Object
{
public:
	DataBuilder();
	virtual int Build() = 0;
};