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

	///Real[0][0]: min x
	///Real[0][1]: max x
	///Real[1][0]: min y
	///Real[1][1]: max y
	///Real[2][0]: min z
	///Real[2][1]: max z
	using Range = Real[3][2];
}
