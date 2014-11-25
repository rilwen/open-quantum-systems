#ifndef __MATH_INTEGR_INTEGRATION_UTILS_H
#define __MATH_INTEGR_INTEGRATION_UTILS_H

namespace rql {
	namespace math {
		namespace integr {
			namespace IntegrationUtils {
				template <class Function> struct RealPart
				{
					RealPart(Function f);
					Function func;
					double operator()(double x)
					{
						return func(x).real();
					}
				};
				template <class F> RealPart<F>::RealPart(F f)
					: func(f)
				{
				}

				template <class Function> struct ImaginaryPart
				{
					ImaginaryPart(Function f);
					Function func;
					double operator()(double x)
					{
						return func(x).imag();
					}
				};
				template <class F> ImaginaryPart<F>::ImaginaryPart(F f)
					: func(f)
				{
				}
			}
		}
	}
}

#endif // __MATH_INTEGR_INTEGRATION_UTILS_H
