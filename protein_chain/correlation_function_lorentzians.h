#ifndef __CORRELATION_FUNCTION_LORENTZIANS_H
#define __CORRELATION_FUNCTION_LORENTZIANS_H

#include "CorrelationFunction.h"
#include "correlation_function_decomposable.h"
#include "complex_analytic_spectral_density.h"
#include "correlation_function_temperature_calculator.h"
#include "core.h"
#include <Eigen/Core>

//! Sum of a/((w-w0)^2 + b2)
class CorrelationFunctionLorentzians: public CorrelationFunctionDecomposable, public virtual ComplexAnalyticSpectralDensity
{
public:
	PROTEIN_CHAIN_API CorrelationFunctionLorentzians(double w0, double hwhm, double height, double kBT=0, size_t pfdOrder=100u);
	PROTEIN_CHAIN_API CorrelationFunctionLorentzians(const Eigen::VectorXd& w0, const Eigen::VectorXd& hwhm, const Eigen::VectorXd& height, double kBT=0, size_t pfdOrder=100u);
	std::complex<double> operator()(double tau) const;
	PROTEIN_CHAIN_API static double height(double scale, double hwhm);
	PROTEIN_CHAIN_API static double scale(double height, double hwhm);
	const Eigen::VectorXd& omegas() const { return m_w0; }
	const Eigen::VectorXd& gammas() const { return m_gamma; }
	const Eigen::VectorXd& scales() const { return m_scale; }
	double kBT() const { return m_kBT; }
	PROTEIN_CHAIN_API double spectral_density(double omega) const;
	PROTEIN_CHAIN_API size_t nbr_exponents() const;
	PROTEIN_CHAIN_API std::complex<double> scale(size_t k) const;
	PROTEIN_CHAIN_API std::complex<double> exponent(size_t k) const;
	std::complex<double> spectral_density(const std::complex<double>& omega) const;
	double spectral_density_derivative(double omega) const;
	size_t nbrPoles() const;
	std::complex<double> pole(size_t idx) const;
	std::complex<double> residue(size_t idx) const;
private:
	int m_nbr_peaks;
	Eigen::VectorXd m_w0;
	Eigen::VectorXd m_gamma;
	Eigen::VectorXd m_scale;	
	double m_kBT;
	CorrelationFunctionTemperatureCalculator m_temp_calc;
};

#endif // __CORRELATION_FUNCTION_LORENTZIANS_H
