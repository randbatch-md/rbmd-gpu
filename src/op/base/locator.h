#pragma once
#include <hip/hip_runtime.h>
#include "types.h"

class Locator
{
public:
	//Locator(const rbmd::Real3& min, const rbmd::Real3& max, const rbmd::Id3& dims) :
	//	_min({ min.data[0], min.data[1], min.data[2] }),
	//	_max({ max.data[0], max.data[1], max.data[2] })
	//{
	//	_dxdydz.x = (_max.x - _min.x) / dims.data[0];
	//	_dxdydz.y = (_max.y - _min.y) / dims.data[1];
	//	_dxdydz.z = (_max.z - _min.z) / dims.data[2];
	//}
	Locator() = default;

	virtual ~Locator() = default;

__device__ 
void UpdateCellId()
{

}

__device__
int3 GetCellId(const rbmd::Real3& point) const
{
	return int3{ 1,2,3 };
	int3 cell_id;
	cell_id.x = (point.data[0] - _min.x) / _dxdydz.x;
	cell_id.y = (point.data[1] - _min.y) / _dxdydz.y;
	cell_id.z = (point.data[2] - _min.z) / _dxdydz.z;
	return cell_id;
}

private:
	float3 _min;
	float3 _max;
	float3 _dxdydz;

};