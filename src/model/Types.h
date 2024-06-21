#include <array>

namespace rbmd
{
#ifdef USE_DOUBLE
	using Real = double;
#else
	using Real = float;
#endif

#ifdef USE_LONGLONG
	using Id = long long;
#else
	using Id = int;
#endif

	using Range = Real[3][2];

	using vec3 = std::array<Real,3>;
}
