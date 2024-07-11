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

		Real3 operator-(const Real3& real) const
		{
			Real3 real3;
			real3.data[0] = this->data[0] - real.data[0];
			real3.data[1] = this->data[1] - real.data[1];
			real3.data[2] = this->data[2] - real.data[2];
			return real3;
		}

		Real3 operator/(const Real3& real) const
		{
			Real3 real3;
			real3.data[0] = this->data[0] / real.data[0];
			real3.data[1] = this->data[1] / real.data[1];
			real3.data[2] = this->data[2] / real.data[2];
			return real3;
		}



		Real data[3];
	};
}
