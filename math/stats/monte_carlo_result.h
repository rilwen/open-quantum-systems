#ifndef __RQL_MATH_STATS_MC_RESULT_H
#define __RQL_MATH_STATS_MC_RESULT_H

#include "../average.h"
#include "../MathCore.h"

namespace rql {
	namespace math {
		namespace stats {
			class MonteCarloResult
			{
			public:				
				RQL_MATH_API_DEBUG MonteCarloResult();
				RQL_MATH_API_DEBUG void update(double x);
				double value() const { return m_average_value.value(); }
				RQL_MATH_API_DEBUG double variance() const;
				RQL_MATH_API_DEBUG double error() const;
			private:
				Average<double> m_average_value;
				Average<double> m_average_squared_value;
			};
		}
	}
}

#endif // __RQL_MATH_STATS_MC_RESULT_H
