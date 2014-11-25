#include "thermal_states.h"
#include <math/linalg/MatrixFunctions.h>
#include <cmath>
#include <stdexcept>

namespace ThermalStates
{
	std::complex<double> thermal_density_matrix(const Eigen::MatrixXcd& H, double kBT, Eigen::MatrixXcd& rho)
	{
		rho = - H/kBT;
		rql::math::linalg::MatrixFunctions::hermitian_matrix_function<double (*)(double)>(rho, exp, rho);
		const std::complex<double> trace_exp = rho.trace();
		rho /= trace_exp;
		return trace_exp;
	}

	double average_mode_occupation(const double omega, const double kBT)
	{
		if (kBT != 0) {
			return 1.0 / (exp(omega/kBT)-1);
		} else {
			if (omega != 0) {
				return 0;
			} else {
				throw std::domain_error("average_mode_occupation: result undefined for omega == kBT == 0");
			}
		}
	}

	double average_mode_energy(double omega, double kBT)
	{
		if (omega != 0) {
			return omega * average_mode_occupation(omega, kBT);
		} else {
			return 0;
		}
	}
}
