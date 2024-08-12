#pragma once
#include <vector>

template<typename T>
struct SOA3DArrayH
{
	std::vector<T> _x;
	std::vector<T> _y;
	std::vector<T> _z;
};
