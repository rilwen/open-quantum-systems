#ifndef __SPECTRAL_DENSITY_H
#define __SPECTRAL_DENSITY_H

#include "core.h"

class SpectralDensity
{
public:
	PROTEIN_CHAIN_API virtual ~SpectralDensity() {}
	PROTEIN_CHAIN_API virtual double spectral_density(double omega) const = 0;
	PROTEIN_CHAIN_API virtual double spectral_density_derivative(double omega) const = 0;
};

#endif // __SPECTRAL_DENSITY_H
