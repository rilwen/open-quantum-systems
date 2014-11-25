#ifndef __MATH_INTERP_INTERPOLATOR_IMPL_FACTORY_H
#define __MATH_INTERP_INTERPOLATOR_IMPL_FACTORY_H

#include <memory>
#include <vector>
#include "../MathCore.h"

namespace rql {
	namespace interp {
		class InterpolatorImpl;

		struct InterpolatorImplFactory
		{
			RQL_MATH_API_DEBUG static std::shared_ptr<InterpolatorImpl> constant(double y);
			RQL_MATH_API_DEBUG static std::shared_ptr<InterpolatorImpl> constant(double y, double lb);
			RQL_MATH_API_DEBUG static std::shared_ptr<InterpolatorImpl> piecewiseConstant(const std::vector<double>& x, const std::vector<double>& y, bool leftInclusive);
			RQL_MATH_API_DEBUG static std::shared_ptr<InterpolatorImpl> piecewiseLinear(const std::vector<double>& x, const std::vector<double>& y);
			RQL_MATH_API_DEBUG static std::shared_ptr<InterpolatorImpl> piecewiseCubic(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& dy);
			RQL_MATH_API_DEBUG static std::shared_ptr<InterpolatorImpl> akima(const std::vector<double>& x, const std::vector<double>& y);
		};
	}
}

#endif // __MATH_INTERP_INTERPOLATOR_IMPL_FACTORY_H
