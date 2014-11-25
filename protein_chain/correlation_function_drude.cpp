#include "correlation_function_drude.h"
#include "partial_function_decomposition_bose.h"
#include "correlation_function_lorentzians.h"
#include <protein_chain/utils.h>
#include <math/MathUtils.h>
#include <stdexcept>
#include <cassert>
#include <cmath>
#include <stdexcept>

CorrelationFunctionDrude::CorrelationFunctionDrude(const double lambda, const double omega_cutoff, const double kBT, const unsigned int order)
	: m_lambda(lambda), m_omega_cutoff(omega_cutoff), m_kBT(kBT)
{
	PartialFunctionDecompositionBose pfd(order);
	m_temp_calc.init(*this, pfd, kBT);
	/*const double eta = 2*lambda/rql::math::PI;
	const std::complex<double> zc(0.0, -omega_cutoff);	
	std::complex<double> tmp = 0.5 + kBT / zc;
	for (unsigned int j = 0; j < order; ++j) {
		const std::complex<double>& xi = pfd.xi(j);
		std::complex<double> z = 2.0 * sqrt(xi) * kBT;
		if (z.imag() == 0) {
			throw std::runtime_error("Real pole detected");
		}
		if (z.imag() > 0) {
			z *= -1;
		}
		m_scale[j] = spectral_density(z) * std::complex<double>(0.0, -2*kBT);
		m_exponent[j] = z * std::complex<double>(0.0, -1.0);
		tmp += kBT * 2 * zc / (zc*zc - 4.0 * xi);
	}
	m_exponent[order] = -omega_cutoff;
	m_scale[order] = tmp * eta * zc;*/
}

std::complex<double> CorrelationFunctionDrude::operator()(double tau) const
{
	return m_temp_calc.evaluate(tau);
	/*const size_t n = nbr_exponents();
	std::complex<double> result = 0.0;
	const double abs_tau = std::abs(tau);
	for (size_t i = 0; i < n; ++i) {
		result += m_scale[i] * exp(abs_tau * m_exponent[i]);
	}
	if (tau > 0) {
		return result;
	} else {
		return conj(result);
	}*/
}

double CorrelationFunctionDrude::spectral_density(double omega) const
{
	return spectral_density(std::complex<double>(omega)).real();
}

boost::shared_ptr<CorrelationFunctionLorentzians> CorrelationFunctionDrude::approximate_by_lorentzian(double lambda, double omega_cutoff, double kBT)
{
	const double hwhm = omega_cutoff;
	const double scale = 2*lambda*kBT;
	const double height = CorrelationFunctionLorentzians::height(scale, hwhm);
	return boost::shared_ptr<CorrelationFunctionLorentzians>(new CorrelationFunctionLorentzians(0.0, hwhm, height, 0.0));
}

std::complex<double> CorrelationFunctionDrude::spectral_density(const std::complex<double>& omega) const
{
	return 2 * m_lambda * omega * m_omega_cutoff / (omega*omega + m_omega_cutoff*m_omega_cutoff) / rql::math::PI;
}

double CorrelationFunctionDrude::spectral_density_derivative(double omega) const
{
	return 2 * m_lambda * m_omega_cutoff / (omega*omega + m_omega_cutoff*m_omega_cutoff) / rql::math::PI
		- 4 * m_lambda * omega * omega * m_omega_cutoff / pow(omega*omega + m_omega_cutoff*m_omega_cutoff, 2) / rql::math::PI;
}

size_t CorrelationFunctionDrude::nbrPoles() const
{
	return 2;
}

std::complex<double> CorrelationFunctionDrude::pole(const size_t idx) const
{
	assert(idx < 2u);
	if (idx == 0) {
		return std::complex<double>(0.0, -m_omega_cutoff);
	} else {
		return std::complex<double>(0.0, m_omega_cutoff);
	}
}

std::complex<double> CorrelationFunctionDrude::residue(size_t idx) const
{
	assert(idx < 2u);
	return m_lambda * m_omega_cutoff / rql::math::PI;
}

//** HIGH-TEMP APPROXIMATION**

//std::complex<double> CorrelationFunctionDrude::operator()(double tau) const
//{
//	assert(m_kBT > 0);
//	// N == 0
//	std::complex<double> result(2*m_lambda*m_kBT,
//		m_only_real ? 0.0 : - m_lambda*m_omega_cutoff);
//	result *= exp(-m_omega_cutoff*std::abs(tau)) * rql::math::PI;
//	if (tau >= 0) {
//		if (tau > 0) {
//			return result;
//		} else {
//			return result.real();
//		}
//	} else {
//		return conj(result);
//	}
//}

