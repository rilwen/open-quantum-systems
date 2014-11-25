#ifndef __QSD_TEST_GAUSSIAN_CORR_FUNC_H
#define __QSD_TEST_GAUSSIAN_CORR_FUNC_H

#include <protein_chain/CorrelationFunction.h>
#include <stdexcept>

class Gaussian: public CorrelationFunction
{
public:
	std::complex<double> operator()(double tau) const
	{
		return exp(-tau*tau);
	}

	/*std::complex<double> atZero(size_t k) const
	{
		if (k % 2 != 0)
			return 0.0;
		else {
			throw std::runtime_error("Not implemented");
		}
	}

	std::complex<double> covariance(size_t k, double d) const
	{
		throw std::runtime_error("Not implemented");
	}*/
};

#endif // __QSD_TEST_GAUSSIAN_CORR_FUNC_H
