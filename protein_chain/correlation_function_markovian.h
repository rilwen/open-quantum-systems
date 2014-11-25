#ifndef __QSD_CORRELATIONFUNCTIONMARKOVIAN_H
#define __QSD_CORRELATIONFUNCTIONMARKOVIAN_H

#include "CorrelationFunction.h"

class CorrelationFunctionMarkovian: public CorrelationFunction
{
public:
	CorrelationFunctionMarkovian(const std::complex<double>& theta);
	std::complex<double> operator()(double tau) const { return 0.0; }
	bool onlyMarkovian() const { return true; }
	/*std::complex<double> atZero(size_t n) const;
	std::complex<double> covariance(size_t k, double d) const;*/
};

#endif // __QSD_CORRELATIONFUNCTIONMARKOVIAN_H
