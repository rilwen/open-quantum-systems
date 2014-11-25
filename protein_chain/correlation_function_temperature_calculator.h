#ifndef __CORRELATION_FUNCTION_TEMPERATURE_CALCULATOR_H
#define __CORRELATION_FUNCTION_TEMPERATURE_CALCULATOR_H

#include <complex>
#include <vector>
#include "core.h"

class ComplexAnalyticSpectralDensity;
class PartialFunctionDecompositionBose;

//! Decomposes a correlation function using PFD for tau > 0
class CorrelationFunctionTemperatureCalculator
{
public:
	CorrelationFunctionTemperatureCalculator();
	PROTEIN_CHAIN_API_DEBUG CorrelationFunctionTemperatureCalculator(const ComplexAnalyticSpectralDensity& J, const PartialFunctionDecompositionBose& pfd, double kBT);
	PROTEIN_CHAIN_API_DEBUG CorrelationFunctionTemperatureCalculator(const ComplexAnalyticSpectralDensity& J, size_t order, double kBT);
	void init(const ComplexAnalyticSpectralDensity& J, const PartialFunctionDecompositionBose& pfd, double kBT);
	void initMatsubara(const ComplexAnalyticSpectralDensity& J, size_t order, double kBT);
	PROTEIN_CHAIN_API_DEBUG std::complex<double> evaluate(double tau) const;
	size_t size() const { return m_size; }
	const std::complex<double>& scale(size_t i) const { return m_scales[i]; }
	const std::complex<double>& exponent(size_t i) const { return m_exponents[i]; }
private:
	bool pole_relevant(const std::complex<double>& pole) const;
	std::complex<double> evaluateMatsubara(const size_t order, const std::complex<double>& z) const;

	size_t m_size;
	std::vector<std::complex<double> > m_scales;
	std::vector<std::complex<double> > m_exponents;
};

#endif // __CORRELATION_FUNCTION_TEMPERATURE_CALCULATOR_H
