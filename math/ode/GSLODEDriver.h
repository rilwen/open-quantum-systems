#ifndef __MATH_GSL_ODE_DRIVER_H
#define __MATH_GSL_ODE_DRIVER_H

#include <gsl/gsl_odeiv2.h>
#include <boost/shared_ptr.hpp>
#include "GSLODESystem.h"
#include "../MathCore.h"

namespace rql {
	namespace math {
		namespace ode {			

			class GSLODEDriver
			{
			public:
				RQL_MATH_API_DEBUG GSLODEDriver(const GSLODEDriver& other);
				RQL_MATH_API_DEBUG int apply(double& t, double ti, double* y) const;
				RQL_MATH_API_DEBUG int applyFixedStep(double& t, const double h, const size_t n, double* y) const;
				RQL_MATH_API_DEBUG size_t dim() const;
				RQL_MATH_API_DEBUG void setHmin(double hMin);
				RQL_MATH_API_DEBUG void setHmax(double hMax);
				RQL_MATH_API_DEBUG void setNmax(size_t nMax);
				RQL_MATH_API_DEBUG static GSLODEDriver y(const boost::shared_ptr<IGSLODESystem>& sys, const gsl_odeiv2_step_type* T, const double hstart, const double epsabs, const double epsrel);
				RQL_MATH_API_DEBUG static GSLODEDriver yp(const boost::shared_ptr<IGSLODESystem>& sys, const gsl_odeiv2_step_type* T, const double hstart, const double epsabs, const double epsrel);
				RQL_MATH_API_DEBUG static GSLODEDriver standard(const boost::shared_ptr<IGSLODESystem>& sys, const gsl_odeiv2_step_type* T, const double hstart, const double epsabs, const double epsrel, const double a_y, const double a_dydt);
				RQL_MATH_API_DEBUG static GSLODEDriver scaled(const boost::shared_ptr<IGSLODESystem>& sys, const gsl_odeiv2_step_type* T, const double hstart, const double epsabs, const double epsrel, const double a_y, const double a_dydt, const double scale_abs[]);
			private:
				GSLODEDriver(gsl_odeiv2_driver* driver, const boost::shared_ptr<IGSLODESystem>& system);
				boost::shared_ptr<IGSLODESystem> m_system;
				boost::shared_ptr<gsl_odeiv2_driver> m_driver;
			};

			inline size_t GSLODEDriver::dim() const
			{
				return m_system->system().dimension;
			}			
		}
	}
}

#endif // __MATH_GSL_ODE_DRIVER_H
