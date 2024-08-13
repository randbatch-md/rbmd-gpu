#pragma once
#include <thrust/device_vector.h>

template<typename T>
struct SOA3DArrayD
{
	thrust::device_vector<T> _x;
	thrust::device_vector<T> _y;
	thrust::device_vector<T> _z;
};
