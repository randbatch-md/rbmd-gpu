#pragma once
#include "object.h"

class MemoryScheduler : public Object
{
public:
	/**
	 * @brief async memeory host to device
	 * @return error code
	*/
	virtual int asyncMemoryH2D() = 0;

	/**
	 * @brief async memory device to host
	 * @return error code
	*/
	virtual int asyncMemoryD2H() = 0;
};