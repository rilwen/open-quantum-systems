#ifndef __MATH_INTERPOLATOR_IMPL_PIECEWISE_CONSTANT_H
#define __MATH_INTERPOLATOR_IMPL_PIECEWISE_CONSTANT_H

#include "InterpolatorImplPiecewise.h"
#include <algorithm>
#include <cassert>

namespace rql {
	namespace interp {

		class InterpolatorImplPiecewiseConstant: public InterpolatorImplPiecewise
		{
		public:
			// x.size() == y.size() + 1
			// if LeftInclusive == true, evaluate(z) == y[i] when x[i] <= z < x[i+1]
			// if LeftInclusive == false, evaluate(z) == y[i] when x[i] < z <= x[i+1]
			// evaluate(z) is always well-defined for z == x[0] or z == x[x.size()-1]
			template <class V1, class V2>
			InterpolatorImplPiecewiseConstant(const V1& x, const V2& y, bool leftInclusive);
			template <class XItBegin, class XItEnd, class YItBegin, class YItEnd>
			InterpolatorImplPiecewiseConstant(XItBegin xBegin, XItEnd xEnd, YItBegin yBegin, YItEnd yEnd, bool leftInclusive);
			RQL_MATH_API_DEBUG InterpolatorImplPiecewiseConstant(const std::vector<double>& x, const std::vector<double>& y, bool leftInclusive);
			RQL_MATH_API_DEBUG virtual std::shared_ptr<InterpolatorImpl> clone() const;
			RQL_MATH_API_DEBUG virtual InterpolatorImpl& operator+=(double x);
			RQL_MATH_API_DEBUG virtual InterpolatorImpl& operator-=(double x);
			RQL_MATH_API_DEBUG virtual InterpolatorImpl& operator*=(double x);
			RQL_MATH_API_DEBUG virtual InterpolatorImpl& operator/=(double x);
		private:
			RQL_MATH_API_DEBUG virtual double evaluate(double x, unsigned int seg_idx) const;
			void setup();
		private:			
			std::vector<double> m_y;
		};

		

		template <class V1, class V2>
		InterpolatorImplPiecewiseConstant::InterpolatorImplPiecewiseConstant(const V1& x, const V2& y, bool leftInclusive)
			: InterpolatorImplPiecewise(x, leftInclusive), m_y(y.size())
		{
			unsigned int i = 0;
			for (std::vector<double>::iterator it = m_y.begin(); it != m_y.end(); ++it)
			{
				*it = y[i];
				++i;
			}
			setup();
		}

		template <class XItBegin, class XItEnd, class YItBegin, class YItEnd>
		InterpolatorImplPiecewiseConstant::InterpolatorImplPiecewiseConstant(XItBegin xBegin, XItEnd xEnd, YItBegin yBegin, YItEnd yEnd, bool leftInclusive)
			: InterpolatorImplPiecewise(xBegin, xEnd, leftInclusive), m_y(0u)
		{
			std::copy(yBegin, yEnd, std::back_inserter(m_y));
			setup();
		}

		inline void InterpolatorImplPiecewiseConstant::setup()
		{
			assert( x().size() == m_y.size() + 1 );
			assert( m_y.size() > 0 );
		}
	}
}

#endif // __MATH_INTERPOLATOR_IMPL_PIECEWISE_CONSTANT_H