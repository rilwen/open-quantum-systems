#ifndef __MATH_INTERPOLATOR_IMPL_CONSTANT_H
#define __MATH_INTERPOLATOR_IMPL_CONSTANT_H

#include "InterpolatorImpl.h"
#include "../MathCore.h"

namespace rql {
	namespace interp {
		class InterpolatorImplConstant: public InterpolatorImpl
		{
		public:
			RQL_MATH_API_DEBUG InterpolatorImplConstant(double constant);
			RQL_MATH_API_DEBUG InterpolatorImplConstant(double constant, double lowerBound);
			RQL_MATH_API_DEBUG virtual double evaluate(double x) const;
			RQL_MATH_API_DEBUG virtual double lowerBound() const;
			RQL_MATH_API_DEBUG virtual double upperBound() const;
			RQL_MATH_API_DEBUG virtual std::shared_ptr<InterpolatorImpl> clone() const;
			RQL_MATH_API_DEBUG virtual InterpolatorImpl& operator+=(double x);
			RQL_MATH_API_DEBUG virtual InterpolatorImpl& operator-=(double x);
			RQL_MATH_API_DEBUG virtual InterpolatorImpl& operator*=(double x);
			RQL_MATH_API_DEBUG virtual InterpolatorImpl& operator/=(double x);
		private:
			double m_const;
			double m_lower_bound;			
		};

		
	}
}

#endif // __MATH_INTERPOLATOR_IMPL_CONSTANT_H