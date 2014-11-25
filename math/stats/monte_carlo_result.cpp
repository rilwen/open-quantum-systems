#include "monte_carlo_result.h"
#include <algorithm>
#include <cassert>
#include <cmath>

namespace rql {
	namespace math {
		namespace stats {
			MonteCarloResult::MonteCarloResult()
				: m_average_value(0), m_average_squared_value(0)
			{}

			void MonteCarloResult::update(double x)
			{
				m_average_value.update(x);
				m_average_squared_value.update(x*x);
			}

			double MonteCarloResult::error() const
			{
				assert(m_average_squared_value.counter() == m_average_value.counter());
				return sqrt(std::max(0.0, (m_average_squared_value.value() - m_average_value.value()*m_average_value.value()) / m_average_value.counter()));
			}
		}
	}
}
