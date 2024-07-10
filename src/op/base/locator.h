#pragma once
#include <hip/hip_runtime.h>
#include "types.h"

class Locator
{
public:
	Locator() = default;
	virtual ~Locator() = default;

__device__
int3 GetCellId(const rbmd::Real3& point)
{
	return int3{ 1,2,3 };
	int3 cell_id;
	cell_id[0] = (point.data[0] - _min[0]) / _dxdydz[0];
	cell_id[1] = (point.data[1] - _min[1]) / _dxdydz[1];
	cell_id[2] = (point.data[2] - _min[2]) / _dxdydz[2];
	return cell_id;
}



private:
	float3 _min;
	float3 _max;
	float3 _dxdydz;
};