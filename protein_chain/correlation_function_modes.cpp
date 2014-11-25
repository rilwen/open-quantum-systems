#include "correlation_function_modes.h"
#include <cmath>
#include <stdexcept>

CorrelationFunctionModes::CorrelationFunctionModes(const Eigen::VectorXcd& g, const Eigen::VectorXd& omega)
	: m_mod_g(g.size()), m_omega(omega)
{
	if (g.size() != omega.size())
		throw std::domain_error("Incompatible argument vectors");
	for (int i = 0; i < g.size(); ++i)
		m_mod_g[i] = pow(std::abs(g[i]), 2);
}

CorrelationFunctionModes::CorrelationFunctionModes()
{
}

std::complex<double> CorrelationFunctionModes::operator()(double tau) const
{
	std::complex<double> sum(0.0, 0.0);
	for (int i = 0; i < m_mod_g.size(); ++i)
		sum += m_mod_g[i] * exp(std::complex<double>(0.0, - m_omega[i]*tau));
	return sum;
}

//std::complex<double> CorrelationFunctionModes::atZero(size_t n) const
//{
//	std::complex<double> sum(0.0, 0.0);
//	for (int i = 0; i < m_mod_g.size(); ++i)
//		sum += m_mod_g[i] * pow(std::complex<double>(0.0, -m_omega[i]), n);
//	return sum;
//}
//
//std::complex<double> CorrelationFunctionModes::covariance(size_t k, double d) const
//{
//	std::complex<double> sum(0.0, 0.0);
//	const double a = (k - 0.5)*d;
//	const double b = (k + 0.5)*d;
//	for (int i = 0; i < m_mod_g.size(); ++i) {
//		sum += m_mod_g[i] * exp(std::complex<double>(0, -m_omega[i]*k*d)) * pow(sin_x_over_x(m_omega[i]*d/2), 2);
//	}
//	return sum * d * d;
//}

