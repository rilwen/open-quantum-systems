#ifndef __CORRELATION_FUNCTION_DECOMPOSABLE_H
#define __CORRELATION_FUNCTION_DECOMPOSABLE_H

#include "CorrelationFunction.h"
#include <complex>

//! Correlation function decomposed for tau > 0 as alpha(tau > 0) = sum_{k=1}^N S[k] * exp(W[k]*tau)
//! alpha(-tau) = conj(tau); if sum_k S[k] has a non-zero imaginary part, it has to be zeroed manually
class CorrelationFunctionDecomposable: public CorrelationFunction
{
public:
	//! Number of terms in the sum
	virtual size_t nbr_exponents() const = 0;
	//! k-th scale S[k]
	virtual std::complex<double> scale(size_t k) const = 0;
	//! k-th exponent W[k]
	virtual std::complex<double> exponent(size_t k) const = 0;
};

#endif // __CORRELATION_FUNCTION_DECOMPOSABLE_H
