#ifndef __MATH_GSL_ODE_SYSTEM_H
#define __MATH_GSL_ODE_SYSTEM_H

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <boost/shared_ptr.hpp>

namespace rql {
	namespace math {
		namespace ode {

			class IGSLODESystem
			{
			public:
				~IGSLODESystem() {}
			private:
				virtual const gsl_odeiv2_system& system() const = 0;
				friend class GSLODEDriver;
			};

			template <class F> class GSLODESystem: public IGSLODESystem
			{
			public:
				static boost::shared_ptr<IGSLODESystem> build(size_t dim, F func);
			private:
				GSLODESystem(size_t dim, F func);
				static int odeFunc(double t, const double* y, double* f, void* params);											
				const gsl_odeiv2_system& system() const { return m_system; }
				gsl_odeiv2_system m_system;
				F m_f;
			};

			template <class F> boost::shared_ptr<IGSLODESystem> GSLODESystem<F>::build(size_t dim, F func)
			{
				return boost::shared_ptr<IGSLODESystem>(new GSLODESystem<F>(dim, func));
			}

			template <class F> GSLODESystem<F>::GSLODESystem(size_t dim, F func)
				: m_f(func)
			{
				m_system.dimension = dim;
				m_system.jacobian = NULL;
				m_system.params = &m_f;
				m_system.function = odeFunc;
			}

			template <class F> int GSLODESystem<F>::odeFunc(double t, const double* y, double* f, void* params)
			{
				F* ode_functor = static_cast<F*>(params);
				return (*ode_functor)(t, y, f);
			}			
		}
	}
}

#endif // __MATH_GSL_ODE_SYSTEM_H
