#include "correlation_function_lorentzians.h"
#include "partial_function_decomposition_bose.h"
#include <stdexcept>
#include <cmath>
#include <complex>
#include <math/MathUtils.h>

static const double SUMMATION_EPSILON = 1E-12;

CorrelationFunctionLorentzians::CorrelationFunctionLorentzians(double w0, double hwhm, double height, double kBT, size_t pfdOrder)
	: m_nbr_peaks(1), m_w0(1), m_gamma(1), m_scale(1), m_kBT(kBT)
{
	if (hwhm <= 0.0) {
		throw std::domain_error("HWHM must be positive");
	}
	m_w0[0] = w0;
	m_gamma[0] = hwhm;
	m_scale[0] = scale(height, hwhm);
	if (kBT > 0) {
		//use Matsubara
		m_temp_calc.initMatsubara(*this, pfdOrder, kBT); 
		/*//use PFD
		PartialFunctionDecompositionBose pfd(pfdOrder);
		m_temp_calc.init(*this, pfd, kBT);*/
	}
}

CorrelationFunctionLorentzians::CorrelationFunctionLorentzians(const Eigen::VectorXd& w0, const Eigen::VectorXd& hwhm, const Eigen::VectorXd& height, double kBT, size_t pfdOrder)
	: m_nbr_peaks(w0.size()), m_w0(w0.size()), m_gamma(w0.size()), m_scale(w0.size()), m_kBT(kBT)
{
	m_w0 = w0;
	if (w0.size() != hwhm.size() || w0.size() != height.size())
		throw std::domain_error("Sizes not consistent");
	for (int i = 0; i < m_nbr_peaks; ++i) {		
		m_gamma[i] = hwhm[i];
		if (hwhm[i] <= 0.0) {
			throw std::domain_error("HWHM must be positive");
		}
		m_scale[i] = scale(height[i], hwhm[i]);
	}
	if (kBT > 0) {
		//use Matsubara
		m_temp_calc.initMatsubara(*this, pfdOrder, kBT); 
		/*//use PFD
		PartialFunctionDecompositionBose pfd(pfdOrder);
		m_temp_calc.init(*this, pfd, kBT);*/
	}
}

std::complex<double> CorrelationFunctionLorentzians::operator()(double tau) const
{		
	if (m_kBT==0) {
		std::complex<double> sum(0);
		// calculate first assuming tau >= 0
		double abs_tau = std::abs(tau);
		for (int i = 0; i < m_nbr_peaks; ++i) {
			sum += m_scale[i]*exp(std::complex<double>( - m_gamma[i]*abs_tau, -m_w0[i]*abs_tau));
		}
		if (tau > 0) {
			return sum;
		} else {
			return conj(sum);
		}
	} else {
		return m_temp_calc.evaluate(tau);
	}
}

double CorrelationFunctionLorentzians::height(double scale, double hwhm)
{
	if (hwhm <= 0)
			throw std::domain_error("HWHM must be positive");
	return scale / hwhm / rql::math::PI;
}

double CorrelationFunctionLorentzians::scale(double height, double hwhm)
{
	if (hwhm <= 0)
			throw std::domain_error("HWHM must be positive");
	return height * hwhm * rql::math::PI;
}

double CorrelationFunctionLorentzians::spectral_density(double omega) const
{
	double ofTheJedi = 0;

	for (int i = 0; i < m_nbr_peaks; ++i) {
		const double gamma = m_gamma[i];
		ofTheJedi += m_scale[i] / rql::math::PI / gamma / (1 + pow((omega - m_w0[i]) / gamma, 2));
	}

	return ofTheJedi;
}

double CorrelationFunctionLorentzians::spectral_density_derivative(double omega) const
{
	double ofTheJedi = 0;

	for (int i = 0; i < m_nbr_peaks; ++i) {
		const double gamma = m_gamma[i];
		const double dw = omega - m_w0[i];
		ofTheJedi -= 2 * m_scale[i] * gamma * dw / pow(dw*dw + gamma*gamma, 2) / rql::math::PI;
	}

	return ofTheJedi;
}

size_t CorrelationFunctionLorentzians::nbr_exponents() const
{
	if (m_kBT != 0) {
		return m_temp_calc.size();
	} else {
		return scales().size();
	}
}

std::complex<double> CorrelationFunctionLorentzians::scale(size_t k) const
{
	if (m_kBT != 0) {
		return m_temp_calc.scale(k);
	} else {
		return scales()[k];
	}
}

std::complex<double> CorrelationFunctionLorentzians::exponent(size_t k) const
{
	if (m_kBT != 0) {
		return m_temp_calc.exponent(k);
	} else {
		return std::complex<double>(-gammas()[k],-omegas()[k]);
	}
}

std::complex<double> CorrelationFunctionLorentzians::spectral_density(const std::complex<double>& omega) const
{
	std::complex<double> ofTheJedi = 0;

	for (int i = 0; i < m_nbr_peaks; ++i) {
		const double gamma = m_gamma[i];
		ofTheJedi += m_scale[i] / rql::math::PI / gamma / (1.0 + pow((omega - m_w0[i]) / gamma, 2));
	}

	return ofTheJedi;
}

size_t CorrelationFunctionLorentzians::nbrPoles() const
{
	return 2*m_nbr_peaks;
}

std::complex<double> CorrelationFunctionLorentzians::pole(size_t idx) const
{
	const size_t peak_idx = idx / 2;
	if (idx % 2 == 0) {
		return std::complex<double>(m_w0[peak_idx], -m_gamma[peak_idx]);
	} else {
		return std::complex<double>(m_w0[peak_idx], m_gamma[peak_idx]);
	}
}

std::complex<double> CorrelationFunctionLorentzians::residue(size_t idx) const
{
	const size_t peak_idx = idx / 2;
	if (idx % 2 == 0) {
		return std::complex<double>(0.0, m_scale[peak_idx] / (2 * rql::math::PI));
	} else {
		return std::complex<double>(0.0, -m_scale[peak_idx] / (2 * rql::math::PI));
	}
}

