#pragma once
#include <cassert>
#include "../MathCore.h"

namespace rql { namespace math { namespace integr {

	/** Maps a univariate integrand from [-infty, infty] to [-1, 1] or [0, infty] to [0, 1].
	 * @tparam F Function object which has operator()(double) defined.
	 * @tparam V Type returned by the function.
	 */
	template <typename F, typename V>
	struct DomainConverterObj
	{
		inline DomainConverterObj(F & f)
		: fun(f)
		{
		}

		/** Calculate the function value. The returned values can be integrated over the finite range to
		 * calculate the corresponding indefinite integral.
		 * @param x Parameter from [-1, 1] or [0, 1] range.
		 * @return Function value for the argument translated to [-infty,infty] or [infty,infty] range times the Jacobian.
		 */
		inline V operator()(double x) const
		{
			assert(x <= 1 && x >= -1);
			const double t = 1 - (x > 0 ? x : -x);
			const double t2 = t*t;
			if (t2 != 0) {
				return fun(x / t) / t2;
			} else {
				return fun(x / t);
			}
		}

		F& fun;
	};

	struct DomainConverter
	{
		// Convert an integration point from the finite [-1, 1] range to [-infty, infty] range.
		RQL_MATH_API_DEBUG static double from_finite_range(double x)
		{
			return x / (1 - (x > 0 ? x : -x));
		}
		// Convert an integration point from the [-infty, infty] range to the finite [-1, 1] range.
		RQL_MATH_API_DEBUG static double to_finite_range(double y);

		RQL_MATH_API_DEBUG static double lower_finite_range() { return - 1; }
		RQL_MATH_API_DEBUG static double upper_finite_range() { return 1; }
	};

}}}