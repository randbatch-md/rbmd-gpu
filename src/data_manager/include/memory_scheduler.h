#pragma once
#include <hip/hip_runtime.h>
#include "object.h"

class MemoryScheduler : public Object
{
public:
	/**
	 * @brief async memeory host to device
	 * @return error code
	*/
	virtual bool asyncMemoryH2D() = 0;

	/**
	 * @brief async memory device to host
	 * @return error code
	*/
	virtual bool asyncMemoryD2H() = 0;
};