#include "correlation_function_markovian.h"
#include <stdexcept>

CorrelationFunctionMarkovian::CorrelationFunctionMarkovian(const std::complex<double>& theta)
	: CorrelationFunction(theta)
{
}

//std::complex<double> CorrelationFunctionMarkovian::atZero(size_t n) const
//{
//	throw std::runtime_error("Not supported");
//}
//
//std::complex<double> CorrelationFunctionMarkovian::covariance(size_t k, double d) const
//{
//	if (k == 0)
//		return d;
//	else
//		return 0;
//}

