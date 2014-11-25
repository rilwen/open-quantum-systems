#ifndef __MATH_GSL_ODE_SOLVER_H
#define __MATH_GSL_ODE_SOLVER_H

#include <gsl/gsl_odeiv2.h>
#include <vector>
#include <Eigen/Core>
#include "../MathCore.h"

namespace rql {
	namespace math {
		namespace ode {

			class GSLODEDriver;

			struct GSLODESolver
			{
				//! first element of x must contain the initial condition
				RQL_MATH_API_DEBUG static void solve(const GSLODEDriver& driver, const std::vector<double>& t, std::vector<Eigen::VectorXd>& x);
			};			
		}
	}
}

#endif // __MATH_GSL_ODE_SOLVER_H
