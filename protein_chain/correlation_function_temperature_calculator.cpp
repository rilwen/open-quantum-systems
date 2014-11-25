#include "correlation_function_temperature_calculator.h"
#include "complex_analytic_spectral_density.h"
#include "partial_function_decomposition_bose.h"
#include <math/MathUtils.h>
#include <boost/math/special_functions/fpclassify.hpp>
#include <stdexcept>
#include <cassert>
#include <cmath>
#include <sstream>
#include <iostream>

CorrelationFunctionTemperatureCalculator::CorrelationFunctionTemperatureCalculator()
	: m_size(0)
{
}

CorrelationFunctionTemperatureCalculator::CorrelationFunctionTemperatureCalculator(const ComplexAnalyticSpectralDensity& J, const PartialFunctionDecompositionBose& pfd, double kBT)
{
	init(J, pfd, kBT);
}

CorrelationFunctionTemperatureCalculator::CorrelationFunctionTemperatureCalculator(const ComplexAnalyticSpectralDensity& J, size_t order, double kBT)
{
	initMatsubara(J, order, kBT);
}

void CorrelationFunctionTemperatureCalculator::init(const ComplexAnalyticSpectralDensity& J, const PartialFunctionDecompositionBose& pfd, double kBT)
{
	if (kBT <= 0.0) {
		throw std::domain_error("Temperature must be positive");
	}
	const double beta = 1.0 / kBT;
	const std::complex<double> residue_multiplier(0.0, -2.0*rql::math::PI); // if tau > 0, use minus sign because the contour is counterclockwise
	static const std::complex<double> MINUS_IMAGINARY(0.0, -1.0);
	std::vector<std::complex<double> > relevant_J_poles;
	std::vector<std::complex<double> > relevant_J_residues;
	relevant_J_poles.reserve(J.nbrPoles());
	relevant_J_residues.reserve(J.nbrPoles());
	for (size_t i = 0; i < J.nbrPoles(); ++i) {
		const std::complex<double> pole_i = J.pole(i);
		if (pole_relevant(pole_i)) {
			relevant_J_poles.push_back(pole_i);
			relevant_J_residues.push_back(J.residue(i));
		}
	}
	m_size = pfd.order() + relevant_J_poles.size();
	m_scales.resize(m_size);
	m_exponents.resize(m_size);
	// do poles of J first
	const size_t nbr_relevant_J_poles = relevant_J_poles.size();
	for (size_t i = 0; i < nbr_relevant_J_poles; ++i) {
		m_scales[i] = residue_multiplier * relevant_J_residues[i] * pfd.evaluate(relevant_J_poles[i]/kBT); // Does not give good results for very small kBT
		m_exponents[i] = MINUS_IMAGINARY * relevant_J_poles[i];
	}
	// skip the pole of PFD in 0.0 so that the correlation function vanishes at infinity
	for (size_t i = 0; i < pfd.order(); ++i) {
		std::complex<double> pole = 2.0*sqrt(pfd.xi(i))*kBT;
		if (!pole_relevant(pole)) {
			pole *= -1;
		}
		const size_t idx = nbr_relevant_J_poles + i;
		m_scales[idx] = residue_multiplier * kBT * J.spectral_density(pole);
		m_exponents[idx] = MINUS_IMAGINARY * pole;
	}

	/*for (size_t i = 0; i < m_size; ++i) {
		std::cout << "Peak " << i << " " << m_scales[i] << " " << m_exponents[i] << "\n";
	}*/

	//// Hack: zero the imaginary parts of all scales
	//for (size_t i = 0; i < m_size; ++i) {
	//	m_scales[i] = m_scales[i].real();
	//}
}

void CorrelationFunctionTemperatureCalculator::initMatsubara(const ComplexAnalyticSpectralDensity& J, const size_t order, double kBT)
{
	if (kBT <= 0.0) {
		throw std::domain_error("Temperature must be positive");
	}
	const double beta = 1.0 / kBT;
	const std::complex<double> residue_multiplier(0.0, -2.0*rql::math::PI); // if tau > 0, use minus sign because the contour is counterclockwise
	static const std::complex<double> MINUS_IMAGINARY(0.0, -1.0);
	std::vector<std::complex<double> > relevant_J_poles;
	std::vector<std::complex<double> > relevant_J_residues;
	relevant_J_poles.reserve(J.nbrPoles());
	relevant_J_residues.reserve(J.nbrPoles());
	for (size_t i = 0; i < J.nbrPoles(); ++i) {
		const std::complex<double> pole_i = J.pole(i);
		if (pole_relevant(pole_i)) {
			relevant_J_poles.push_back(pole_i);
			relevant_J_residues.push_back(J.residue(i));
		}
	}
	m_size = order + relevant_J_poles.size();
	m_scales.resize(m_size);
	m_exponents.resize(m_size);
	// do poles of J first
	const size_t nbr_relevant_J_poles = relevant_J_poles.size();
	for (size_t i = 0; i < nbr_relevant_J_poles; ++i) {
		m_scales[i] = residue_multiplier * relevant_J_residues[i] * evaluateMatsubara(order, relevant_J_poles[i]/kBT); // Does not give good results for very small kBT
		m_exponents[i] = MINUS_IMAGINARY * relevant_J_poles[i];
	}
	// skip the pole of Matsubara decomposition in 0.0 so that the correlation function vanishes at infinity
	for (size_t i = 0; i < order; ++i) {
		std::complex<double> pole(0.0, 2.0*rql::math::PI*(i+1)*kBT);
		if (!pole_relevant(pole)) {
			pole *= -1;
		}
		const size_t idx = nbr_relevant_J_poles + i;
		m_scales[idx] = residue_multiplier * kBT * J.spectral_density(pole);
		m_exponents[idx] = MINUS_IMAGINARY * pole;
	}

	/*for (size_t i = 0; i < m_size; ++i) {
		std::cout << "Peak " << i << " " << m_scales[i] << " " << m_exponents[i] << "\n";
	}*/
}

bool CorrelationFunctionTemperatureCalculator::pole_relevant(const std::complex<double>& pole) const
{
	if (pole.imag() == 0.0) {
		std::stringstream ss;
		ss << "Cannot handle entirely real poles: " << pole;
		throw std::runtime_error(ss.str().c_str());
	}
	// integral will be times exp(-i*omega*tau) for tau > 0, hence we integrate over the Im(z) < 0 semi-circle
	return pole.imag() < 0.0;
}

std::complex<double> CorrelationFunctionTemperatureCalculator::evaluate(double tau) const
{
	const double used_tau = std::abs(tau);
	std::complex<double> sum(0.0);
	for (size_t i = 0; i < m_size; ++i) {
		sum += m_scales[i] * exp(m_exponents[i] * used_tau);
	}
	if (!(boost::math::isfinite(sum.real()) && boost::math::isfinite(sum.imag()))) {
		std::stringstream ss;
		ss << "Bad result: " << sum;
		throw std::runtime_error(ss.str().c_str());
	}
	//return sum.real();
	if (tau == 0.0) {
		return sum.real();
	} else {
		if (tau > 0) {
			return sum;
		} else {
			return conj(sum);
		}
	}
}

std::complex<double> CorrelationFunctionTemperatureCalculator::evaluateMatsubara(const size_t order, const std::complex<double>& z) const
{
	std::complex<double> sum = 0.5 + 1.0/z;
	for (size_t i = 0; i < order; ++i)
		sum += 2.0*z/(z*z + 4.0*rql::math::PI*rql::math::PI*(i+1)*(i+1));

	return sum;
}
