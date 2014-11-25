#ifndef __CORRELATION_FUNCTION_DRUDE_H
#define __CORRELATION_FUNCTION_DRUDE_H

#include "correlation_function_decomposable.h"
#include "complex_analytic_spectral_density.h"
#include "correlation_function_temperature_calculator.h"
#include "core.h"
#include <Eigen/Core>
#include <boost/shared_ptr.hpp>

class CorrelationFunctionLorentzians;

//! J(w) = 2*lambda*w*w_C / (w*w + w_C*w_C) / PI
class CorrelationFunctionDrude: public CorrelationFunctionDecomposable, public virtual ComplexAnalyticSpectralDensity
{
public:
	PROTEIN_CHAIN_API CorrelationFunctionDrude(double lambda, double omega_cutoff, double kBT, unsigned int order);
	PROTEIN_CHAIN_API std::complex<double> operator()(double tau) const;	
	PROTEIN_CHAIN_API double spectral_density(double omega) const;
	size_t nbr_exponents() const { return m_temp_calc.size(); };
	std::complex<double> scale(size_t k) const { return m_temp_calc.scale(k); }
	std::complex<double> exponent(size_t k) const { return m_temp_calc.exponent(k); }
	//! for high kBT
	PROTEIN_CHAIN_API static boost::shared_ptr<CorrelationFunctionLorentzians> approximate_by_lorentzian(double lambda, double omega_cutoff, double kBT);
	std::complex<double> spectral_density(const std::complex<double>& omega) const;
	double spectral_density_derivative(double omega) const;
	size_t nbrPoles() const;
	std::complex<double> pole(size_t idx) const;
	std::complex<double> residue(size_t idx) const;
private:
	double m_lambda;
	double m_omega_cutoff;
	double m_kBT;
	CorrelationFunctionTemperatureCalculator m_temp_calc;
};

#endif // __CORRELATION_FUNCTION_DRUDE_H
