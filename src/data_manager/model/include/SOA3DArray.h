#pragma once
#include <vector>

template<typename T>
struct SOA3DArray
{
	std::vector<T> _x;
	std::vector<T> _y;
	std::vector<T> _z;
};
