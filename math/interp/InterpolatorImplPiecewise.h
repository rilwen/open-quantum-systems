#ifndef __MATH_INTERP_INTERPOLATOR_IMPL_PIECEWISE_H
#define __MATH_INTERP_INTERPOLATOR_IMPL_PIECEWISE_H

#include <algorithm>
#include <cassert>
#include <vector>
#include "InterpolatorImpl.h"
#include "../MathCore.h"

namespace rql {
	namespace interp {

		// Type of interpolator which is defined piecewise.
		class InterpolatorImplPiecewise: public InterpolatorImpl
		{
		public:
			// if LeftInclusive == true, find such i that x[i] <= z < x[i+1]
			// if LeftInclusive == false, find such i that x[i] < z <= x[i+1]
			// evaluate(z) is always well-defined for z == x[0] or z == x[x.size()-1]
			template <class V>
			InterpolatorImplPiecewise(const V& x, bool leftInclusive);
			template <class ItBegin, class ItEnd>
			InterpolatorImplPiecewise(ItBegin xBegin, ItEnd xEnd, bool leftInclusive);
			RQL_MATH_API_DEBUG InterpolatorImplPiecewise(const std::vector<double>& x, bool leftInclusive);
			RQL_MATH_API_DEBUG virtual double evaluate(double x) const;
			RQL_MATH_API_DEBUG virtual double lowerBound() const;
			RQL_MATH_API_DEBUG virtual double upperBound() const;
		private:
			// x - initial argument
			// seg_idx - index of the LEFT segment boundary in m_x
			RQL_MATH_API_DEBUG virtual double evaluate(double x, unsigned int seg_idx) const = 0;
			void setup();
		protected:
			const std::vector<double>& x() const;
			bool left_inclusive() const;
			unsigned int nbr_segments() const;
		private:
			std::vector<double> m_x;
			const bool m_left_inclusive;
			unsigned int m_nbr_segments;
		};

		template <class V>
		InterpolatorImplPiecewise::InterpolatorImplPiecewise(const V& x, bool leftInclusive)
			: m_x(x.size()), m_left_inclusive(leftInclusive)
		{
			unsigned int i = 0;
			for (std::vector<double>::iterator it = m_x.begin(); it != m_x.end(); ++it)
			{
				*it = x[i];
				++i;
			}
			setup();
		}

		template <class ItBegin, class ItEnd>
		InterpolatorImplPiecewise::InterpolatorImplPiecewise(ItBegin xBegin, ItEnd xEnd, bool leftInclusive)
			: m_x(0u), m_left_inclusive(leftInclusive)
		{
			std::copy(xBegin, xEnd, std::back_inserter(m_x));
			setup();
		}

		inline void InterpolatorImplPiecewise::setup()
		{
			if (m_x.size() <= 1)
				throw std::runtime_error("m_x size too small");
			m_nbr_segments = m_x.size() - 1;
		}

		inline const std::vector<double>& InterpolatorImplPiecewise::x() const
		{
			return m_x;
		}

		inline unsigned int InterpolatorImplPiecewise::nbr_segments() const
		{
			return m_nbr_segments;
		}
	}
}

#endif // __MATH_INTERP_INTERPOLATOR_IMPL_PIECEWISE_H