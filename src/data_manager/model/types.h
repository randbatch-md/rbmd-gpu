#pragma once

namespace rbmd
{
#if USE_DOUBLE
	using Real = double;
#else
	using Real = float;
#endif

#if USE_64BIT_IDS
	using Id = long long;
#else
	using Id = int;
#endif

	using Range = Real[3][2];
}
