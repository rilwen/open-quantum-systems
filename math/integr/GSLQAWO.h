#ifndef __MATH_INTEGR_GSLQAWO_H
#define __MATH_INTEGR_GSLQAWO_H

#include <gsl/gsl_integration.h>
#include "../MathCore.h"

namespace rql {
	namespace math {
		namespace integr {
			class GSLQAWO
			{
			public:
				RQL_MATH_API_DEBUG GSLQAWO(double omega, double L, gsl_integration_qawo_enum type, size_t n);
				RQL_MATH_API_DEBUG ~GSLQAWO();
				RQL_MATH_API_DEBUG int setOmega(double omega);
				RQL_MATH_API_DEBUG int setL(double L);
				RQL_MATH_API_DEBUG int setType(gsl_integration_qawo_enum type);
				RQL_MATH_API_DEBUG int set(double omega, double L, gsl_integration_qawo_enum type);
				RQL_MATH_API_DEBUG void integrateImpl(gsl_function* f, double a, double epsabs, double epsrel, size_t limit, double& result, double& abserr);			
				template <class Function> void integrate(Function integrand, double a, double epsabs, double epsrel, size_t limit, double& result, double& abserr);			
			private:
				template <class Function> static double integrand(double x, void* params);			
			private:
				double m_omega;
				double m_L;
				gsl_integration_qawo_enum m_type;
				size_t m_n;
				gsl_integration_workspace* m_wksp;
				gsl_integration_qawo_table* m_table;
			};

			template <class F> void GSLQAWO::integrate(F fun, double a, double epsabs, double epsrel, size_t limit, double& result, double& abserr)
			{
				gsl_function f;
				f.function = &integrand<F>;
				f.params = &fun;
				integrateImpl(&f, a, epsabs, epsrel, limit, result, abserr);
			}

			template <class F> double GSLQAWO::integrand(double x, void* params)
			{
				F* functor = static_cast<F*>(params);
				return (*functor)(x);
			}
		}
	}
}

#endif // __MATH_INTEGR_GSLQAWO_H
