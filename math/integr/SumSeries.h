#ifndef __RQL_MATH_INTEGR_SUM_SERIES_H
#define __RQL_MATH_INTEGR_SUM_SERIES_H

#include <math/Accumulator.h>
#include <cmath>

namespace rql { namespace math { namespace integr {
	template <class V>
	class SumSeries
	{
	public:		
		SumSeries(double epsilon = 1E-15);
		template <class S>
		V sum(S series, const V& initial_value, int start, int step) const;
		template <class S>
		V sum(S series, const V& initial_value, unsigned int nbr_initial_steps, int start, int step) const;
	private:
		double m_epsilon;
	};

	template <class V>
	SumSeries<V>::SumSeries(double epsilon)
		: m_epsilon(epsilon)
	{
	}

	template <class V>
	template <class S>
	V SumSeries<V>::sum(S series, const V& initial_value, int start, int step) const
	{
		Accumulator<V> acc;
		acc = initial_value;
		int idx = start;
		V old_sum(0);
		while (idx >= start) {
			acc += series(idx);
			const V new_sum = acc.sum();
			if (fabs(new_sum - old_sum) < (fabs(new_sum) + 1) * m_epsilon) {
				break;
			}
			old_sum = new_sum;
			idx += step;
		}
		return acc.sum();
	}

	template <class V>
	template <class S>
	V SumSeries<V>::sum(S series, const V& initial_value, unsigned int nbr_initial_steps, int start, int step) const
	{
		Accumulator<V> acc;
		acc = initial_value;
		for (size_t i = 0; i < nbr_initial_steps; ++i) {
			acc += series(start + i*step);
		}
		acc += sum(series, start + step*nbr_initial_steps, step);
		return acc.sum();
	}

}}}

#endif // __RQL_MATH_INTEGR_SUM_SERIES_H
