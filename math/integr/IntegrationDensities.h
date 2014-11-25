#ifndef __MATH_INTEGR_INTEGRATION_DENSITIES_H
#define __MATH_INTEGR_INTEGRATION_DENSITIES_H

#include "../MathCore.h"

namespace rql {
	namespace integr {
		namespace IntegrationDensities {
			//! w(x) = a
			class Constant
			{
			public:
				RQL_MATH_API_DEBUG Constant(double a) { m_a = a; }
				RQL_MATH_API_DEBUG double operator()(double x0, double x1, unsigned int k) const;
			private:
				double m_a;
			};

			//! w(x) = sin(omega*x)
			class Sine
			{
			public:
				RQL_MATH_API_DEBUG Sine(double omega) { m_omega = omega; }
				RQL_MATH_API_DEBUG double operator()(double x0, double x1, unsigned int k) const;
			private:
				double m_omega;
			};

			// f(x) = cos(omega*x)
			class Cosine
			{
			public:
				RQL_MATH_API_DEBUG Cosine(double omega) { m_omega = omega; }
				RQL_MATH_API_DEBUG double operator()(double x0, double x1, unsigned int k) const;
			private:
				double m_omega;
			};
		}		
	}
}

#endif // __MATH_INTEGR_INTEGRATION_DENSITIES_H
