#include "FourierIntegral.h"
#include "IntegrationDensities.h"
#include <stdexcept>
#include <cassert>

namespace rql {
	namespace integr {
		FourierIntegral::FourierIntegral(double omega, const std::vector<double>& x)
			: m_size(x.size()), m_wi_sin(x, IntegrationDensities::Sine(omega)), m_wi_cos(x, IntegrationDensities::Cosine(omega))
			, m_y_r(x.size()), m_y_i(x.size())
		{
			if (x.size() < 2)
				throw std::domain_error("At least two X points needed");
		}

		FourierIntegral::FourierIntegral(double omega, double x0, double x1, size_t nbrSegments)
			: m_size(nbrSegments + 1), m_wi_sin(WeightedIntegral::gen_xs(x0, x1, nbrSegments), IntegrationDensities::Sine(omega)), m_wi_cos(WeightedIntegral::gen_xs(x0, x1, nbrSegments), IntegrationDensities::Cosine(omega))
			, m_y_r(nbrSegments + 1), m_y_i(nbrSegments + 1)
		{
		}

		std::complex<double> FourierIntegral::integrate(const std::vector<std::complex<double> >& y)
		{
			assert(m_size == y.size());
			std::vector<double>::iterator rit = m_y_r.begin();
			std::vector<double>::iterator iit = m_y_i.begin();
			for (std::vector<std::complex<double> >::const_iterator yit = y.begin(); yit != y.end(); ++yit,++rit,++iit) {
				assert(rit != m_y_r.end());
				assert(iit != m_y_i.end());
				*rit = yit->real();
				*iit = yit->imag();
			}
			return integrate(m_y_r, m_y_i);			
		}

		std::complex<double> FourierIntegral::integrate(const std::vector<double>& yr, const std::vector<double>& yi)
		{
			const double icR = m_wi_cos.integrate(yr);
			const double icI = m_wi_cos.integrate(yi);
			const double isR = m_wi_sin.integrate(yr);
			const double isI = m_wi_sin.integrate(yi);
			return std::complex<double>( icR - isI, icI + isR );
		}
	}
}
