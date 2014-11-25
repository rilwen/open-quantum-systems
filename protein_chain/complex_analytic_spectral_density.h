#ifndef __COMPLEX_ANALYTIC_SPECTRAL_DENSITY_H
#define __COMPLEX_ANALYTIC_SPECTRAL_DENSITY_H

#include "spectral_density.h"
#include <complex>

class ComplexAnalyticSpectralDensity: public SpectralDensity
{
public:
	virtual std::complex<double> spectral_density(const std::complex<double>& omega) const = 0;
	virtual size_t nbrPoles() const = 0;
	virtual std::complex<double> pole(size_t idx) const = 0;
	virtual std::complex<double> residue(size_t idx) const = 0;
};

#endif // __COMPLEX_ANALYTIC_SPECTRAL_DENSITY_H
