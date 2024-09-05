#pragma once
#include <string>
#include "object.h"
class MDData;

class InBuilder : public Object
{
public:
	InBuilder();
	int Build() override;

protected:
	MDData& _md_data;
};