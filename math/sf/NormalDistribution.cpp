#include <cassert>
#include <limits>
#include "NormalDistribution.h"

namespace rql
{
	namespace math
	{
		namespace sf
		{

			const double calerf_xneg = -26.628;
			const double calerf_xinf = std::numeric_limits<double>::max();

			static inline double calerf_return(double x, bool calculate_erfc, double result)
			{
				if (calculate_erfc) {
					return x < 0 ? 2 - result : result;
				} else {
					result = (0.5 - result) + 0.5;
					return (x < 0) ? - result : result;
				}
			}

			const double calerf_thresh = 0.46875;


			// Ported from http://www.netlib.org/specfun/erf
			// Fortran code packet by W. J. Cody
			// Removed support for error function * exp(-x^2)

			double calerf(const double x, const bool calculate_erfc)
			{
				static const double invsqrtpi = 0.5641895835477562869480794516;

				// Portable
				static const double xsmall = std::numeric_limits<double>::epsilon();

				// Valid for double precision in IEEE
				static const double xbig = 26.543;

				static const double a[] = {3.16112374387056560E0
					, 1.13864154151050156E2
					, 3.77485237685302021E2
					, 3.20937758913846947E3
					, 1.85777706184603153E-1};

				static const double b[] = {2.36012909523441209E1, 2.44024637934444173E2, 1.28261652607737228E3, 2.84423683343917062E3};

				static const double c[] = {5.64188496988670089E-1, 8.88314979438837594E0,
					6.61191906371416295E1, 2.98635138197400131E2,
					8.81952221241769090E2, 1.71204761263407058E3,
					2.05107837782607147E3, 1.23033935479799725E3,
					2.15311535474403846E-8};

				static const double d[] = {1.57449261107098347E1, 1.17693950891312499E2,
					5.37181101862009858E2, 1.62138957456669019E3,
					3.29079923573345963E3, 4.36261909014324716E3,
					3.43936767414372164E3, 1.23033935480374942E3};

				static const double p[] = {3.05326634961232344E-1, 3.60344899949804439E-1,
					1.25781726111229246E-1, 1.60837851487422766E-2,
					6.58749161529837803E-4, 1.63153871373020978E-2};

				static const double q[] = {2.56852019228982242E0, 1.87295284992346047E0,
					5.27905102951428412E-1, 6.05183413124413191E-2,
					2.33520497626869185E-3};

				const double y = std::abs(x);
				register double xnum;
				register double xden;

				// Evaluate  erf  for  |X| <= 0.46875
				if (y <= calerf_thresh) {
					const double ysq = y > xsmall ? y*y : 0;
					xnum = (a[4] * ysq + a[0]) * ysq;
					xden = (ysq + b[0]) * ysq;
					xnum = (xnum + a[1]) * ysq;
					xden = (xden + b[1]) * ysq;
					xnum = (xnum + a[2]) * ysq;
					xden = (xden + b[2]) * ysq;
					return calculate_erfc ? 1 - x * (xnum + a[3]) / (xden + b[3]) : x * (xnum + a[3]) / (xden + b[3]);
				}		// Evaluate  erfc  for 0.46875 <= |X| <= 4.0
				else if (y <= 4) {
					xnum = (c[8] * y + c[0]) * y;
					xden = (y + d[0]) * y;
					xnum = (xnum + c[1]) * y;
					xden = (xden + d[1]) * y;
					xnum = (xnum + c[2]) * y;
					xden = (xden + d[2]) * y;
					xnum = (xnum + c[3]) * y;
					xden = (xden + d[3]) * y;
					xnum = (xnum + c[4]) * y;
					xden = (xden + d[4]) * y;
					xnum = (xnum + c[5]) * y;
					xden = (xden + d[5]) * y;
					xnum = (xnum + c[6]) * y;
					xden = (xden + d[6]) * y;
					return calerf_return(x, calculate_erfc, ((xnum + c[7]) / (xden + d[7])) * exp(-y*y));
				}		// Evaluate  erfc  for |X| > 4.0
				else {
					if (y >= xbig) {
						return calerf_return(x, calculate_erfc, 0.0);
					}
					register const double ysq = 1 / (y * y);
					xnum = (p[5]*ysq + p[0]) * ysq;
					xden = (ysq + q[0]) * ysq;
					xnum = (xnum + p[1]) * ysq;
					xden = (xden + q[1]) * ysq;
					xnum = (xnum + p[2]) * ysq;
					xden = (xden + q[2]) * ysq;
					xnum = (xnum + p[3]) * ysq;
					xden = (xden + q[3]) * ysq;
					return calerf_return(x, calculate_erfc, ((invsqrtpi - (ysq * (xnum + p[4]) / (xden + q[4]))) / y) * exp(-y*y));
				}

			}

			static const double p_low = 0.02425;
			static const double p_high = 1 - p_low;

			// Implementation of the algorithm described in http://home.online.no/~pjacklam/notes/invnorm/
			double normsinv(double p)
			{
				assert( p >= 0 );
				assert( p <= 1 );
				int sign = 1;
				if (p > 0.5) {
					p = 1 - p;
					sign = -1;
				}
				if (p == 0) {
					return -sign * std::numeric_limits<double>::infinity();
				}

				static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02,
					-2.759285104469687e+02, 1.383577518672690e+02,
					-3.066479806614716e+01, 2.506628277459239e+00};
				static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02,
					-1.556989798598866e+02, 6.680131188771972e+01,
					-1.328068155288572e+01};
				static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01,
					-2.400758277161838e+00, -2.549732539343734e+00,
					4.374664141464968e+00, 2.938163982698783e+00};
				static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01,
					2.445134137142996e+00, 3.754408661907416e+00};
				double x;
				if (p < p_low) {
					const double q = sqrt(-2*log(p));
					x = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
						((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
				} else {
					const double q = p - 0.5;
					const double r = q*q;
					x = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
						(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
				}
				const double e = normcdf(x) - p;
				if (x > -37) {
					const double u = e * 2.506628274631000502415765285 * exp(0.5*x*x);
					return sign * (x - u/(1 + 0.5*x*u));
				} else {
					return x;
				}
			}
		}
	}
}
