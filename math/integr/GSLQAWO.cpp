#include "GSLQAWO.h"
#include <stdexcept>
#include <gsl/gsl_errno.h>
#include "../GSLError.h"

namespace rql {
	namespace math {
		namespace integr {
			GSLQAWO::GSLQAWO(double omega, double L, gsl_integration_qawo_enum type, size_t n)
				: m_omega(omega), m_L(L), m_type(type), m_n(n)
			{
				m_wksp = gsl_integration_workspace_alloc(n);
				if (!m_wksp)
					throw std::bad_alloc();
				m_table = gsl_integration_qawo_table_alloc(omega, L, type, n);
				if (!m_table) {
					gsl_integration_workspace_free(m_wksp);
					throw std::bad_alloc();
				}
			}

			GSLQAWO::~GSLQAWO()
			{
				gsl_integration_workspace_free(m_wksp);
				gsl_integration_qawo_table_free(m_table);
			}

			void GSLQAWO::integrateImpl(gsl_function* f, double a, double epsabs, double epsrel, size_t limit, double& result, double& abserr)
			{
				const int err = gsl_integration_qawo(f, a, epsabs, epsrel, limit, m_wksp, m_table, &result, &abserr);
				if (err != GSL_SUCCESS)
					throw GSLError("Error in GSL QAWO", err);
			}

			int GSLQAWO::setOmega(double omega)
			{
				return set(omega, m_L, m_type);
			}

			int GSLQAWO::setL(double L)
			{
				return gsl_integration_qawo_table_set_length(m_table, L);
			}

			int GSLQAWO::setType(gsl_integration_qawo_enum type)
			{
				return set(m_omega, m_L, type);
			}

			int GSLQAWO::set(double omega, double L, gsl_integration_qawo_enum type)
			{
				return gsl_integration_qawo_table_set(m_table, omega, L, type);
			}
		}
	}
}
