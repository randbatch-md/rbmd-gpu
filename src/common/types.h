#pragma once

namespace rbmd {
#if USE_DOUBLE
using Real = double;
#else
using Real = float;
#endif

#if USE_64BIT_IDS
using Id = unsigned long;
using Int = long long;
#else
using Id = int;
// using Int = int;
#endif
}  // namespace rbmd
