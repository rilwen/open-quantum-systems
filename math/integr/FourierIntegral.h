#ifndef __MATH_INTEGR_FOURIER_INTEGRAL_H
#define __MATH_INTEGR_FOURIER_INTEGRAL_H

#include "WeightedIntegral.h"
#include "../MathCore.h"
#include <complex>
#include <vector>

namespace rql {
	namespace integr {
		class FourierIntegral
		{
		public:
			RQL_MATH_API_DEBUG FourierIntegral(double omega, const std::vector<double>& x);
			RQL_MATH_API_DEBUG FourierIntegral(double omega, double x0, double x1, size_t nbrSegments);
			RQL_MATH_API_DEBUG std::complex<double> integrate(const std::vector<std::complex<double> >& y);
			RQL_MATH_API_DEBUG std::complex<double> integrate(const std::vector<double>& yr, const std::vector<double>& yi);
			//! @param[in] integrand Function to integrate with given density; returns complex numbers
			template <class F> std::complex<double> integrate(F integrand);
		private:
			FourierIntegral& operator=(const FourierIntegral&); // not implemented
			const size_t m_size;
			WeightedIntegral m_wi_sin;
			WeightedIntegral m_wi_cos;
			std::vector<double> m_y_r;
			std::vector<double> m_y_i;
		};

		template <class F> std::complex<double> FourierIntegral::integrate(F integrand)
		{
			std::vector<double>::iterator rit = m_y_r.begin();
			std::vector<double>::iterator iit = m_y_i.begin();
			for (std::vector<double>::const_iterator xit = m_wi_cos.x().begin(); xit != m_wi_cos.x().end(); ++xit,++rit,++iit) {
				assert(rit != m_y_r.end());
				assert(iit != m_y_i.end());
				const std::complex<double> y = integrand(*xit);
				*rit = y.real();
				*iit = y.imag();
			}
			return integrate(m_y_r, m_y_i);			
		}
	}
}

#endif // __MATH_INTEGR_FOURIER_INTEGRAL_H
