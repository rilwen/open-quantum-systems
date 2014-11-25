#pragma once

#include <cmath>
#include "../MathCore.h"

namespace rql
{
	namespace math
	{
		namespace sf
		{
			// Don't call this directly unless you know what you're doing.
			RQL_MATH_API double calerf(const double x, const bool calculate_erfc);

			// ERF function; most precise around 0
			inline double erf(double x)
			{
				return calerf(x, false);
			}

			// ERFC function; most precise in the tails.
			inline double erfc(double x)
			{
				return calerf(x, true);
			}

			// Normal density function.
			inline double normpdf(double x)
			{
				return 0.398942280401433 * exp(-0.5*x*x);
			}

			inline double normpdf(double x, double mean, double sigma)
			{
				return normpdf((x - mean)/sigma)/sigma;
			}

			// Normal Gaussian CDF
			inline double normcdf(const double x)
			{
				// This ensures max precision.
				return 0.5*erfc(-0.7071067811865475244008443621 * x);
			}

			// Inverse normal Gaussian CDF
			RQL_MATH_API double normsinv(double p);
		}
	}
}