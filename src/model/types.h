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

		Real3 operator-(const Real3& id3)
		{
			Real3 real3;
			real3[0] = this->data[0] - id3[0];
			real3[1] = this->data[1] - id3[1];
			real3[2] = this->data[2] - id3[2];
			return real3;
		}

		Id3 operator/(const Real3& id3)
		{
			Id3 id3;
			id3.data[0] = this->data[0] / id3[0];
			id3.data[1] = this->data[1] / id3[1];
			id3.data[2] = this->data[2] / id3[2];
			return id3;
		}

		Real data[3];
	};
}
