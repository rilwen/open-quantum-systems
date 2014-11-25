#include "CorrelationFunction.h"
#include <stdexcept>
#include <cassert>

CorrelationFunction::CorrelationFunction()
	: m_has_markovian(false), m_markovian(0.0)
{
}

CorrelationFunction::CorrelationFunction(const std::complex<double>& markovian)
	: m_has_markovian(true), m_markovian(markovian)
{
}

std::complex<double> CorrelationFunction::average_exponential(const std::complex<double>& w, double a, double b)
{
	assert( b > a );
	const std::complex<double> x1 = w*a;
	const std::complex<double> x2 = w*b;
	const std::complex<double> dx = x2 - x1;
	if (std::abs(dx) > 1E-8)
		return (exp(x2) - exp(x1))/dx;
	else
		return exp(0.5*x1 + 0.5*x2);
}
