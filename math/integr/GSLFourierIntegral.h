#ifndef __MATH_INTEGR_GSL_FOURIER_INTEGRAL_H
#define __MATH_INTEGR_GSL_FOURIER_INTEGRAL_H

#include "GSLQAWO.h"
#include "IntegrationUtils.h"
#include <complex>
#include <cassert>

namespace rql {
	namespace math {
		namespace integr {
			class GSLFourierIntegral
			{
			public:
				RQL_MATH_API_DEBUG GSLFourierIntegral(double omega, double L, size_t n);
				template <class Function> void integrate(Function integrand, double a, double epsabs, double epsrel, size_t limit, std::complex<double>& result, std::complex<double>& abserr);
			private:
				GSLQAWO m_sin;
				GSLQAWO m_cos;
			};

			template <class F> void GSLFourierIntegral::integrate(F integrand, double a, double epsabs, double epsrel, size_t limit, std::complex<double>& result, std::complex<double>& abserr)
			{
				IntegrationUtils::RealPart<F> rp(integrand);
				IntegrationUtils::ImaginaryPart<F> ip(integrand);
				double sr,cr,si,ci;
				double err_sr, err_cr, err_si, err_ci;
				abserr = 0;
				m_sin.integrate(rp, a, epsabs/2, epsrel, limit, sr, err_sr);
				assert(err_sr >= 0);
				m_cos.integrate(rp, a, epsabs/2, epsrel, limit, cr, err_cr);
				assert(err_cr >= 0);
				m_sin.integrate(ip, a, epsabs/2, epsrel, limit, si, err_si);
				assert(err_si >= 0);
				m_cos.integrate(ip, a, epsabs/2, epsrel, limit, ci, err_ci);
				assert(err_ci >= 0);
				result = std::complex<double>(cr - si, ci + sr);
				abserr = std::complex<double>(err_cr + err_si, err_ci + err_sr);
			}
		}
	}
}

#endif // __MATH_INTEGR_GSL_FOURIER_INTEGRAL_H
