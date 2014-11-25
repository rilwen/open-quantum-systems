#ifndef __ALPHA_COVARIANCE_CALCULATOR_H
#define __ALPHA_COVARIANCE_CALCULATOR_H

#include <complex>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <Eigen/Core>
#include "core.h"

class SingleAlphaCovarianceCalculator;
class CorrelationFunction;

class AlphaCovarianceCalculator
{
public:
	QSD_API_DEBUG AlphaCovarianceCalculator(const std::vector<boost::shared_ptr<const CorrelationFunction> >& alpha, size_t nbrTimePts, double d);
	QSD_API_DEBUG AlphaCovarianceCalculator(boost::shared_ptr<const CorrelationFunction> alpha, const size_t dim, const size_t nbrTimePts, double d);
	double delta() const { return m_d; }
	size_t dim() const { return m_alpha_cov.rows(); }
	size_t nbr_time_pts() const { return m_alpha_cov.cols(); }
	//! return cov(z_m(k*d)^*, z_m(l*d))
	//! m < dim()
	//! k,l < nbrTimePts()
	QSD_API_DEBUG std::complex<double> covariance(size_t m, size_t k, size_t l) const;
	//! return cov(Real(z_m(k*d), Real(z_m(l*d))) == cov(Imag(z_m(k*d), Imag(z_m(l*d)))
	//! m < dim()
	//! k,l < nbrTimePts()
	QSD_API_DEBUG double covariance_RR_II(size_t m, size_t k, size_t l) const;
	//! return cov(Real(z_m(k*d), Imag(z_m(l*d))) == - cov(Imag(z_m(k*d), Real(z_m(l*d)))
	//! m < dim()
	//! k,l < nbrTimePts()
	QSD_API_DEBUG double covariance_RI(size_t m, size_t k, size_t l) const;
	//! Return the calculator for m-th alpha function
	//! m < dim()
	QSD_API_DEBUG const SingleAlphaCovarianceCalculator single_alpha_calculator(size_t m) const;
private:
	double m_d;
	Eigen::MatrixXcd m_alpha_cov;
};

class SingleAlphaCovarianceCalculator
{
public:
	double delta() const { return m_acc.delta(); }
	size_t nbr_time_pts() const { return m_acc.nbr_time_pts(); }
	//! return cov(z_m(k*d)^*, z_m(l*d))
	//! k,l < nbrTimePts()
	std::complex<double> covariance(size_t k, size_t l) const { return m_acc.covariance(m_m, k, l); }
	//! return cov(Real(z(k*d), Real(z(l*d))) == cov(Imag(z(k*d), Imag(z(l*d)))
	//! k,l < nbrTimePts()
	double covariance_RR_II(size_t k, size_t l) const { return m_acc.covariance_RR_II(m_m, k, l); }
	//! return cov(Real(z(k*d), Imag(z(l*d))) == - cov(Imag(z(k*d), Real(z(l*d)))
	//! k,l < nbrTimePts()
	double covariance_RI(size_t k, size_t l) const { return m_acc.covariance_RI(m_m, k, l); }
private:
	SingleAlphaCovarianceCalculator(const AlphaCovarianceCalculator& acc, size_t m);
	const AlphaCovarianceCalculator& m_acc;
	const size_t m_m;
	friend class AlphaCovarianceCalculator;
};


#endif // __ALPHA_COVARIANCE_CALCULATOR_H
