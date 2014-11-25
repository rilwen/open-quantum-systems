#include "absorption_spectrum_calculator_time_dependent.h"
#include <cmath>
#include <stdexcept>
#include <math/MathUtils.h>

AbsorptionSpectrumCalculatorTimeDependent::AbsorptionSpectrumCalculatorTimeDependent(double w0, double w1, double dw)
	: m_w0(w0), m_dw(dw), m_N(static_cast<size_t>(floor((w1 - w0)/dw) + 1))
{
}

void AbsorptionSpectrumCalculatorTimeDependent::calculateAbsorptionSpectrum(const std::vector<Eigen::VectorXcd>& states, const double dt, std::vector<double>& spectrum) const
{
	if (dt <= 0)
		throw std::domain_error("Dt must be positive");	
	std::vector<std::complex<double> > dots(states.size());
	std::vector<double> dot_times(states.size());
	dots[0] = 1;
	dot_times[0] = 0;
	for (size_t j = 1; j < states.size(); ++j) {
		dots[j] = states.front().dot(states[j]);
		dot_times[j] = j*dt;
	}
	calculateAbsorptionSpectrum(dot_times, dots, spectrum);
}

void AbsorptionSpectrumCalculatorTimeDependent::calculateAbsorptionSpectrum(const std::vector<double>& scalProdTimes, const std::vector<std::complex<double> >& scalProds, std::vector<double>& spectrum) const
{
	if (scalProdTimes.size() < 2)
		throw std::domain_error("Too few data points");
	if (scalProdTimes.size() != scalProds.size())
		throw std::domain_error("Wrong sizes");
	spectrum.resize(m_N);
	const size_t m = scalProdTimes.size();
	for (size_t i = 0; i < m_N; ++i) {
		const double w = m_w0 + i*m_dw;
		double sum = 0;
		for (size_t j = 0; j < m; ++j) {
			const double integrand = (exp(std::complex<double>(0.0, w*scalProdTimes[j])) * scalProds[j]).real();
			double weight = 0;
			if (j > 0)
				weight += 0.5*(scalProdTimes[j] - scalProdTimes[j-1]);
			if (j < m - 1)
				weight += 0.5*(scalProdTimes[j+1] - scalProdTimes[j]);
			if (weight <= 0)
				throw std::runtime_error("Weights must be positive");
			sum += weight * integrand;
		}
		spectrum[i] = sum * 4 * rql::math::PI;
	}
}

