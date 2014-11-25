#include "WeightedIntegral.h"
#include "../interp/AkimaApproximation.h"
#include "../interp/PolynomialInterpolation.h"
#include <boost/array.hpp>
#include <stdexcept>

namespace rql {
	namespace integr {
		std::vector<double> WeightedIntegral::gen_xs(double x0, double x1, const size_t nbrSegments)
		{
			if (x1 <= x0)
				throw std::domain_error("Bad [x0,x1] range");
			if (!nbrSegments)
				throw std::domain_error("At least 1 segment needed");
			std::vector<double> x(nbrSegments + 1);
			const double dx = (x1 - x0) / nbrSegments;
			x[0] = x0;
			for (size_t i = 1; i < nbrSegments; ++i)
				x[i] = x0 + i*dx;
			x[nbrSegments] = x1;
			return x;
		}

		WeightedIntegral::WeightedIntegral(const std::vector<double>& x, const std::vector<boost::array<double,ORDER> >& weights)
			: m_size(x.size()), m_x(x), m_weights(weights), m_y(x.size()), m_dy(x.size())
		{
			if (m_size != weights.size() + 1)
				throw std::domain_error("X and WEIGHTS size mismatch");
		}

		double WeightedIntegral::integrate(const std::vector<double>& y)
		{
			assert(m_size == y.size());
			interp::AkimaApproximation<double>::calculate(m_x, y, m_dy);
			boost::array<double,2u> y0;
			boost::array<double,2u> yx;
			boost::array<double,4u> a;
			yx[0] = y[0];
			yx[1] = m_dy[0];
			double integral = 0;
			for (size_t i = 1; i < m_size; ++i) {
				const double x = m_x[i] - m_x[i - 1];
				y0 = yx;
				yx[0] = y[i];
				yx[1] = m_dy[i];
				interp::PolynomialInterpolation::interpolate(x, y0, yx, a);
				const boost::array<double,ORDER>& w = m_weights[i - 1];
				integral += a[0]*w[0] + a[1]*w[1] + a[2]*w[2] + a[3]*w[3];
			}
			return integral;
		}
	}
}
