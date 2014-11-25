#include "GSLODESolver.h"
#include <cassert>
#include <gsl/gsl_errno.h>
#include <sstream>
#include "GSLODEDriver.h"

namespace rql {
	namespace math {
		namespace ode {

			void GSLODESolver::solve(const GSLODEDriver& driver, const std::vector<double>& t, std::vector<Eigen::VectorXd>& x)
			{
				assert(t.size() == x.size());
				assert(t.size() > 1);
				assert(x[0].size() == driver.dim());
				for (size_t i = 1; i < x.size(); ++i) {
					double prev_t = t[i - 1];
					x[i] = x[i - 1];
					const int status = driver.apply(prev_t, t[i], x[i].data());
					if (status != GSL_SUCCESS) {
						std::stringstream ss;
						ss << "GSL error " << status << ": " << gsl_strerror(status);
						throw std::runtime_error(ss.str().c_str());
					}
				}
			}			

		}
	}
}

