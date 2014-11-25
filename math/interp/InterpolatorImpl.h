#ifndef __MATH_INTERPOLATOR_IMPL_H
#define __MATH_INTERPOLATOR_IMPL_H

#include <memory>
#include <boost/utility.hpp>
#include "../MathCore.h"

namespace rql {
	namespace interp {

		class InterpolatorImpl: boost::noncopyable
		{
		public:
			RQL_MATH_API_DEBUG virtual ~InterpolatorImpl();
			RQL_MATH_API_DEBUG virtual double evaluate(double x) const = 0;
			RQL_MATH_API_DEBUG virtual double lowerBound() const = 0;
			RQL_MATH_API_DEBUG virtual double upperBound() const = 0;

			//! Perform a deep copy.
			RQL_MATH_API_DEBUG virtual std::shared_ptr<InterpolatorImpl> clone() const = 0;

			//! Add x to interpolated function.
			RQL_MATH_API_DEBUG virtual InterpolatorImpl& operator+=(double x) = 0;

			//! Subtract x from interpolated function.
			RQL_MATH_API_DEBUG virtual InterpolatorImpl& operator-=(double x) = 0;

			//! Multiply the interpolated function by x.
			RQL_MATH_API_DEBUG virtual InterpolatorImpl& operator*=(double x) = 0;

			//! Divide the interpolated function by x.
			RQL_MATH_API_DEBUG virtual InterpolatorImpl& operator/=(double x) = 0;
		};
	}
}

#endif // __MATH_INTERPOLATOR_IMPL_H