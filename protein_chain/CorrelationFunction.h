#ifndef __CORRELATION_FUNCTION_H
#define __CORRELATION_FUNCTION_H
#include <Eigen/Core>
#include <cmath>
#include "core.h"

//! alpha(tau) function
class CorrelationFunction
{
public:
	PROTEIN_CHAIN_API_DEBUG CorrelationFunction();
	CorrelationFunction(const std::complex<double>& markovian);
	virtual ~CorrelationFunction() {}
	virtual std::complex<double> operator()(double tau) const = 0;
	bool hasMarkovian() const { return m_has_markovian; }
	const std::complex<double>& markovianScale() const { return m_markovian; }
	virtual bool onlyMarkovian() const { return false; }
	////! n-th right derivative at zero
	//virtual std::complex<double> atZero(size_t n) const = 0;
	////! int_{-d/2}^{d/2} dt int_{-d/2}^{d/2} ds alpha(k*d + t - s)
	//virtual std::complex<double> covariance(size_t k, double d) const = 0;
protected:
	//! Calculate int_a^b exp(i*x*w) dx / (b - a)
	//! a < b
	static std::complex<double> average_exponential(const std::complex<double>& w, double a, double b);
	static double sin_x_over_x(double x) { return (x != 0) ? sin(x)/x : 1; }
private:
	bool m_has_markovian;
	std::complex<double> m_markovian;
};

#endif // __CORRELATION_FUNCTION_H
