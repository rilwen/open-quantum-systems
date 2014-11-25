#ifndef __MATH_INTERP_POLYNOMIAL_INTERPOLATION_H
#define __MATH_INTERP_POLYNOMIAL_INTERPOLATION_H

#include <cassert>
#include <vector>

namespace rql {
	namespace interp {
		/** Calculates polynomial interpolation on a range [0, x] */
		namespace PolynomialInterpolation
		{
			//! @param x Must be != 0
			//! @param y0 Value of interpolated function in 0
			//! @param yx Value of interpolated function in x
			//! @param a0 0-th order coefficient
			//! @param a1 1-st order coefficient
			void interpolate(double x, const double y0, const double yx, double& a0, double& a1);

			template <template <class X,size_t> class Array, size_t N>
			struct do_interpolate
			{
				static void apply(double x, const Array<double,N> y0, const Array<double,N> yx, Array<double,static_cast<size_t>(2*N)>& a);
			};

			//! @param x Must be != 0
			//! @param y0 Moments of interpolated function in 0
			//! @param yx Moments of interpolated function in x
			//! @param a Interpolating polynomial coefficients
			//! @tparam N Number of function moments (value, 1st derivative, 2nd derivative, ...)
			template <template <class X,size_t> class Array, size_t N>
			void interpolate(double x, const Array<double,N> y0, const Array<double,N> yx, Array<double,static_cast<size_t>(2*N)>& a)
			{
				do_interpolate<Array,N>::apply(x, y0, yx, a);
			}

			inline void interpolate(double x, const double y0, const double yx, double& a0, double& a1)
			{
				assert( x != 0 );
				a0 = y0;
				a1 = (yx - y0) / x;
			}

			template <template <class X,size_t> class Array>
			struct do_interpolate<Array,1u>
			{
				static void apply(double x, const Array<double,1u> y0, const Array<double,1u> yx, Array<double,2u>& a)
				{
					interpolate(x, y0[0], yx[0], a[0], a[1]);
				}
			};
			

			template <template <class X,size_t> class Array>
			struct do_interpolate<Array,2u>
			{
				static void apply(double x, const Array<double,2u> y0, const Array<double,2u> yx, Array<double,4u>& a)
				{
					assert( x != 0 );
					a[0] = y0[0];
					a[1] = y0[1];
					a[2] = (-yx[1] - 2*y0[1] + 3*(yx[0] - y0[0])/x)/x;
					a[3] = (yx[1] + y0[1] + 2*(y0[0] - yx[0])/x)/x/x;
				}
			};
			
		};		
	}
}

#endif // __MATH_INTERP_POLYNOMIAL_INTERPOLATION_H

