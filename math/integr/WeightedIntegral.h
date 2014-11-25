#ifndef __MATH_INTEGR_WEIGHTED_INTEGRAL_H
#define __MATH_INTEGR_WEIGHTED_INTEGRAL_H

#include <algorithm>
#include <stdexcept>
#include <vector>
#include <boost/array.hpp>
#include "../MathCore.h"

namespace rql {
	namespace integr {
		//! @tparam ValueType Function value type
		class WeightedIntegral
		{
		public:
			RQL_MATH_API_DEBUG static const size_t ORDER = 4;
			RQL_MATH_API_DEBUG static std::vector<double> gen_xs(double x0, double x1, size_t nbrSegments);
			template <class Density> WeightedIntegral(const std::vector<double>& x, Density density);
			RQL_MATH_API_DEBUG WeightedIntegral(const std::vector<double>& x, const std::vector<boost::array<double,ORDER> >& weights);
			RQL_MATH_API_DEBUG double integrate(const std::vector<double>& y);
			//! @param[in] integrand Function to integrate with given density
			template <class F> double integrate(F integrand);
			const std::vector<double>& x() const { return m_x; }
		private:
			WeightedIntegral& operator=(const WeightedIntegral&); // not implemented
			const size_t m_size;
			std::vector<double> m_x;
			std::vector<boost::array<double,ORDER> > m_weights;
			std::vector<double> m_y;
			std::vector<double> m_dy;
		};

		template <class Density> WeightedIntegral::WeightedIntegral(const std::vector<double>& x, Density density)
			: m_size(x.size()), m_x(x), m_y(x.size()), m_dy(x.size())
		{
			if (m_size < 2)
				throw std::domain_error("X size too low");
			m_weights.resize(m_size - 1);
			for (size_t i = 1; i < m_size; ++i) {
				for (size_t k = 0; k < ORDER; ++k) {
					m_weights[i - 1][k] = density(x[i - 1], x[i], k);
				}
			}
		}

		template <class F>
		double WeightedIntegral::integrate(F integrand)
		{
			std::transform(m_x.begin(), m_x.end(), m_y.begin(), integrand);
			return integrate(m_y);
		}
	}
}

#endif // __MATH_INTEGR_WEIGHTED_INTEGRAL_H
