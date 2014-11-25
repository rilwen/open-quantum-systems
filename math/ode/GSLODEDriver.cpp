#include "GSLODEDriver.h"
#include <cassert>
#include "GSLODESystem.h"
#include <stdexcept>

namespace rql {
	namespace math {
		namespace ode {

			GSLODEDriver::GSLODEDriver(const GSLODEDriver& other)
				: m_driver(other.m_driver)
			{
			}

			GSLODEDriver::GSLODEDriver(gsl_odeiv2_driver* driver, const boost::shared_ptr<IGSLODESystem>& system)
				: m_driver(driver, gsl_odeiv2_driver_free), m_system(system)
			{
				m_driver->sys = &m_system->system(); // point to the local copy
			}

			int GSLODEDriver::apply(double& t, double ti, double* y) const
			{
				return gsl_odeiv2_driver_apply(m_driver.get(), &t, ti, y);
			}			

			GSLODEDriver GSLODEDriver::y(const boost::shared_ptr<IGSLODESystem>& sys, const gsl_odeiv2_step_type* T, const double hstart, const double epsabs, const double epsrel)
			{
				return GSLODEDriver(gsl_odeiv2_driver_alloc_y_new(&sys->system(), T, hstart, epsabs, epsrel), sys);
			}

			GSLODEDriver GSLODEDriver::yp(const boost::shared_ptr<IGSLODESystem>& sys, const gsl_odeiv2_step_type* T, const double hstart, const double epsabs, const double epsrel)
			{
				return GSLODEDriver(gsl_odeiv2_driver_alloc_yp_new(&sys->system(), T, hstart, epsabs, epsrel), sys);
			}

			GSLODEDriver GSLODEDriver::standard(const boost::shared_ptr<IGSLODESystem>& sys, const gsl_odeiv2_step_type* T, const double hstart, const double epsabs, const double epsrel, const double a_y, const double a_dydt)
			{
				return GSLODEDriver(gsl_odeiv2_driver_alloc_standard_new(&sys->system(), T, hstart, epsabs, epsrel, a_y, a_dydt), sys);
			}

			GSLODEDriver GSLODEDriver::scaled(const boost::shared_ptr<IGSLODESystem>& sys, const gsl_odeiv2_step_type* T, const double hstart, const double epsabs, const double epsrel, const double a_y, const double a_dydt, const double* scale_abs)
			{
				return GSLODEDriver(gsl_odeiv2_driver_alloc_scaled_new(&sys->system(), T, hstart, epsabs, epsrel, a_y, a_dydt, scale_abs), sys);
			}

			void GSLODEDriver::setHmax(double hMax)
			{
				const int status = gsl_odeiv2_driver_set_hmax(m_driver.get(), hMax);
				if (status != GSL_SUCCESS)
					throw std::runtime_error(gsl_strerror(status));
			}

			void GSLODEDriver::setHmin(double hMin)
			{
				const int status = gsl_odeiv2_driver_set_hmin(m_driver.get(), hMin);
				if (status != GSL_SUCCESS)
					throw std::runtime_error(gsl_strerror(status));
			}

			void GSLODEDriver::setNmax(size_t nMax)
			{
				const int status = gsl_odeiv2_driver_set_nmax(m_driver.get(), nMax);
				if (status != GSL_SUCCESS)
					throw std::runtime_error(gsl_strerror(status));
			}

			int GSLODEDriver::applyFixedStep(double& t, const double h, const size_t n, double* y) const
			{
				return gsl_odeiv2_driver_apply_fixed_step(m_driver.get(), &t, h, n, y);
			}
		}
	}
}
