#pragma once

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

	struct Id3
	{
		Id3() = default;
		Id3(const Id& x, const Id& y, const Id& z)
		{
			data[0] = x;
			data[1] = y;
			data[2] = z;
		}


		Id data[3];
	};

	struct Real3
	{
		Real3() = default;
		Real3(const Real& x, const Real& y, const Real& z)
		{
			data[0] = x;
			data[1] = y;
			data[2] = z;
		}

		Real3 operator-(const Real3& real1, const Real3& real2) const
		{
			Real3 real3;
			real3.data[0] = real1.data[0] - rea2.data[0];
			real3.data[1] = real1.data[1] - rea2.data[1];
			real3.data[2] = real1.data[2] - rea2.data[2];
			return real3;
		}

		Real3 operator/(const Real3& real1, const Real3& real2) const
		{
			Real3 real3;
			real3.data[0] = real1.data[0] / rea2.data[0];
			real3.data[1] = real1.data[1] / rea2.data[1];
			real3.data[2] = real1.data[2] / rea2.data[2];
			return real3;
		}



		Real data[3];
	};
}
