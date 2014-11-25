#include "Brent.h"

namespace rql { namespace math { namespace solvers {
#ifdef NDEBUG
	double root_finder_result(double root_position, double value_at_root, double x1, double f1, double x2, double /*f2*/, double& lower_bound, double& higher_bound)
#else
	double root_finder_result(double root_position, double value_at_root, double x1, double f1, double x2, double f2, double& lower_bound, double& higher_bound)
#endif
	{
		if (value_at_root == 0) {
			lower_bound = root_position;
			higher_bound = root_position;
		} else {
			assert( opposite_signs(f1, f2) );
			assert( f1 != 0 );
			assert( f2 != 0 );
			if (opposite_signs(value_at_root, f1)) {
				lower_bound = std::min(x1, root_position);
				higher_bound = std::max(x1, root_position);
			} else {
				lower_bound = std::min(x2, root_position);
				higher_bound = std::max(x2, root_position);
			}
		}
		return root_position;
	}
}}}
