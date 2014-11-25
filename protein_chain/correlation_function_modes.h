#ifndef __QSD_CORRELATION_FUNCTION_MODES_H
#define __QSD_CORRELATION_FUNCTION_MODES_H

#include "CorrelationFunction.h"
#include <Eigen/Core>

class CorrelationFunctionModes: public CorrelationFunction
{
public:
	CorrelationFunctionModes(const Eigen::VectorXcd& g, const Eigen::VectorXd& omega);
	//! Creates alpha === 0
	CorrelationFunctionModes();
	std::complex<double> operator()(double tau) const;
	bool onlyMarkovian() const { return m_mod_g.size()==0; }
	//std::complex<double> atZero(size_t n) const;
	//std::complex<double> covariance(size_t k, double d) const;
private:
	Eigen::VectorXd m_mod_g;
	Eigen::VectorXd m_omega;
};

#endif // __QSD_CORRELATION_FUNCTION_MODES_H
