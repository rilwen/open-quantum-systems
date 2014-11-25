#include <algorithm>
#include <cassert>
#include <cmath>
#include <cerrno>
#include <limits>
#include "MathUtils.h"

namespace rql { namespace math {

	long gcd(long u, long v)
	{
		if (u == 0 || v == 0) {
			return 0;
		}

		// Both u and v are not zero.

		if (u == v) {
			return u;
		}

		// u != v

		int sign = 1;
		if (u < 0) {
			sign = -sign;
			u = -u;
		}
		if (v < 0) {
			sign = -sign;
			v = -v;
		}

		// u > 0 and v > 0

		if (u > v) {
			std::swap(u, v);
		}

		// u < v

		// Let's go!

		do {
			v = v % u;
			std::swap(u, v);
		}	while (u > 0);

		return sign*v;
	}

	const double PI = 3.14159265358979323846264338327950288;
}}