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
	cell_id[0] = (point.data[0] - _min.x) / _dxdydz.x;
	cell_id[1] = (point.data[1] - _min.y) / _dxdydz.y;
	cell_id[2] = (point.data[2] - _min.z) / _dxdydz.z;
	return cell_id;
}



private:
	float3 _min;
	float3 _max;
	float3 _dxdydz;
};