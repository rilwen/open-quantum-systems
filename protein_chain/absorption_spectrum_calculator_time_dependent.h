#ifndef __TIMEDEPENDENTABSORPTIONSPECTRUMCALCULATOR_H
#define __TIMEDEPENDENTABSORPTIONSPECTRUMCALCULATOR_H

#include <vector>
#include <Eigen/Core>
#include "core.h"

class AbsorptionSpectrumCalculatorTimeDependent
{
public:
	PROTEIN_CHAIN_API AbsorptionSpectrumCalculatorTimeDependent(double w0, double w1, double dw);
	PROTEIN_CHAIN_API void calculateAbsorptionSpectrum(const std::vector<Eigen::VectorXcd>& states, const double dt, std::vector<double>& spectrum) const;
	PROTEIN_CHAIN_API void calculateAbsorptionSpectrum(const std::vector<double>& scalProdTimes, const std::vector<std::complex<double> >& scalProds, std::vector<double>& spectrum) const;
private:
	double m_w0;
	double m_dw;
	size_t m_N;
};

#endif // __TIMEDEPENDENTABSORPTIONSPECTRUMCALCULATOR_H
