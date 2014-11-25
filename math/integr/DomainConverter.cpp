#include "DomainConverter.h"
#include <limits>

namespace rql { namespace math { namespace integr {

	double DomainConverter::to_finite_range(double y)
	{
		if (y == std::numeric_limits<double>::infinity()) {
			return 1;
		} else if (y == -std::numeric_limits<double>::infinity()) {
			return -1;
		} else {
			return y / (1 + (y > 0 ? y : -y));
		}
	}

}}}